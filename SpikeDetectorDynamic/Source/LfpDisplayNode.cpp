/*
    ------------------------------------------------------------------

    This file is part of the Open Ephys GUI
    Copyright (C) 2016 Open Ephys

    ------------------------------------------------------------------

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "LfpDisplayNode.h"
#include "LfpDisplayCanvas.h"
#include <stdio.h>

using namespace LfpViewer;


LfpDisplayNode::LfpDisplayNode()
    : GenericProcessor  ("LFP Spike Viewer")
    , displayGain       (1)
    , bufferLength      (10.0f) //in terms of seconds
    , abstractFifo      (100)
{
    setProcessorType (PROCESSOR_TYPE_SINK);

    displayBuffer = new AudioSampleBuffer (8, 100);

    const int heapSize = 5000;
    arrayOfOnes = new float[heapSize];
    for (int n = 0; n < heapSize; ++n)
    {
        arrayOfOnes[n] = 1;
    }

	subprocessorToDraw = 0;
	numSubprocessors = -1;
}


LfpDisplayNode::~LfpDisplayNode()
{
    delete[] arrayOfOnes;
}


AudioProcessorEditor* LfpDisplayNode::createEditor()
{
    editor = new LfpDisplayEditor (this, true);
    return editor;
}


void LfpDisplayNode::updateSettings()
{

    std::cout << "Setting num inputs on LfpDisplayNode to " << getNumInputs() << std::endl;

	numChannelsInSubprocessor.clear();
    subprocessorSampleRate.clear();

	for (int i = 0; i < getNumInputs(); i++)
	{
        uint32 channelSubprocessor = getDataSubprocId(i);

        numChannelsInSubprocessor.insert({ channelSubprocessor, 0 }); // (if not already there)
        numChannelsInSubprocessor[channelSubprocessor]++;

        subprocessorSampleRate.insert({ channelSubprocessor, getDataChannel(i)->getSampleRate() });
	}
    
    numSubprocessors = numChannelsInSubprocessor.size();

    if (numChannelsInSubprocessor.find(subprocessorToDraw) == numChannelsInSubprocessor.end())
    {
        // subprocessor to draw does not exist
        if (numSubprocessors == 0)
        {
            subprocessorToDraw = 0;
        }
        else
        {
            // there are channels, but none on the current subprocessorToDraw
            // default to the first subprocessor
            subprocessorToDraw = getDataSubprocId(0);
        }
    }

    int numChans = getNumSubprocessorChannels();
    int srate = getSubprocessorSampleRate(subprocessorToDraw);

	std::cout << "Re-setting num inputs on LfpDisplayNode to " << numChans << std::endl;
    if (numChans > 0)
    {
        std::cout << "Sample rate = " << srate << std::endl;
    }

    eventSourceNodes.clear();
    ttlState.clear();

	for (int i = 0; i < eventChannelArray.size(); ++i)
	{
		uint32 sourceId = getEventSourceId(eventChannelArray[i]);
 
		if (!eventSourceNodes.contains(sourceId))
		{
			eventSourceNodes.add(sourceId);
		}
	}

    for (int i = 0; i < eventSourceNodes.size(); ++i)
    {
		std::cout << "Adding channel " << numChans + i << " for event source node " << eventSourceNodes[i] << std::endl;

        ttlState[eventSourceNodes[i]] = 0;
    }

    resizeBuffer();
    
    // update the editor's subprocessor selection display and sample rate
	LfpDisplayEditor * ed = (LfpDisplayEditor*)getEditor();
	ed->updateSubprocessorSelectorOptions();

    //update spike channel
    electrodes.clear();
    for (int i = 0; i < spikeChannelArray.size(); ++i)
    {
        std::cout << "Adding electrode " << std::endl;

        Electrode* elec = new Electrode();
        elec->numChannels = spikeChannelArray[i]->getNumChannels();
        elec->bitVolts = spikeChannelArray[i]->getChannelBitVolts(0); //lets assume all channels have the same bitvolts
        elec->name = spikeChannelArray[i]->getName();
        elec->currentSpikeIndex = 0;
        //elec->mostRecentSpikes.ensureStorageAllocated(displayBufferSize);

        for (int j = 0; j < elec->numChannels; ++j)
        {
            elec->displayThresholds.add(0);
            elec->detectorThresholds.add(0);
        }

        electrodes.add(elec);

    }
}

uint32 LfpDisplayNode::getEventSourceId(const EventChannel* event)
{
    return getProcessorFullId(event->getTimestampOriginProcessor(), event->getTimestampOriginSubProcessor());
}

uint32 LfpDisplayNode::getChannelSourceId(const InfoObjectCommon* chan)
{
    return getProcessorFullId(chan->getSourceNodeID(), chan->getSubProcessorIdx());
}

uint32 LfpDisplayNode::getDataSubprocId(int chan) const
{
    if (chan < 0 || chan >= getTotalDataChannels())
    {
        return 0;
    }

    return getChannelSourceId(getDataChannel(chan));
}

OwnedArray<Electrode>* LfpViewer::LfpDisplayNode::getElectrodes()
{
    return &electrodes;
}



void LfpDisplayNode::setSubprocessor(uint32 sp)
{

	subprocessorToDraw = sp;
    resizeBuffer();
	std::cout << "LfpDisplayNode setting subprocessor to " << sp << std::endl;	
}

uint32 LfpDisplayNode::getSubprocessor() const
{
    return subprocessorToDraw;
}

int LfpDisplayNode::getNumSubprocessorChannels()
{
    if (subprocessorToDraw != 0)
    {
        return numChannelsInSubprocessor[subprocessorToDraw];
    }
    return 0;
}

float LfpDisplayNode::getSubprocessorSampleRate(uint32 subprocId)
{
    auto entry = subprocessorSampleRate.find(subprocId);
    if (entry != subprocessorSampleRate.end())
    {
        return entry->second;
    }
    return 0.0f;
}

bool LfpDisplayNode::resizeBuffer()
{
	int nSamples = (int)getSubprocessorSampleRate(subprocessorToDraw) * bufferLength;
	int nInputs = getNumSubprocessorChannels();

	std::cout << "Resizing buffer. Samples: " << nSamples << ", Inputs: " << nInputs << std::endl;

	if (nSamples > 0 && nInputs > 0)
	{
		abstractFifo.setTotalSize(nSamples);
		displayBuffer->setSize(nInputs + 1, nSamples); // add extra channel for TTLs
		displayBuffer->clear();

		displayBufferIndex.clear();
		displayBufferIndex.insertMultiple(0, 0, nInputs + 1);

		return true;
	}
	else
	{
		return false;
	}
	
}


bool LfpDisplayNode::enable()
{

	if (resizeBuffer())
	{
		LfpDisplayEditor* editor = (LfpDisplayEditor*)getEditor();
		editor->enable();
		return true;
	}
	else
	{
		return false;
	}

}


bool LfpDisplayNode::disable()
{
    LfpDisplayEditor* editor = (LfpDisplayEditor*) getEditor();
    editor->disable();
    return true;
}


void LfpDisplayNode::setParameter (int parameterIndex, float newValue)
{
    editor->updateParameterButtons (parameterIndex);
    //
    //Sets Parameter in parameters array for processor
    parameters[parameterIndex]->setValue (newValue, currentChannel);

    //std::cout << "Saving Parameter from " << currentChannel << ", channel ";

    LfpDisplayEditor* ed = (LfpDisplayEditor*) getEditor();
    if (ed->canvas != 0)
        ed->canvas->setParameter (parameterIndex, newValue);
}


void LfpDisplayNode::handleEvent(const EventChannel* eventInfo, const MidiMessage& event, int samplePosition)
{
    if (Event::getEventType(event) == EventChannel::TTL)
    {
        TTLEventPtr ttl = TTLEvent::deserializeFromMessage(event, eventInfo);
        
        //int eventNodeId = *(dataptr+1);
        const int eventId = ttl->getState() ? 1 : 0;
        const int eventChannel = ttl->getChannel();
        const int eventTime = samplePosition;

        // find sample rate of event channel
        uint32 eventSourceNodeId = getEventSourceId(eventInfo);
        float eventSampleRate = getSubprocessorSampleRate(eventSourceNodeId);

        if (eventSampleRate == 0)
        {
            // shouldn't happen for any real event channel at this point
            return;
        }
        
		//std::cout << "Received event on channel " << eventChannel << std::endl;
		//std::cout << "Copying to channel " << channelForEventSource[eventSourceNodeId] << std::endl;
        
        if (eventId == 1)
        {
            ttlState[eventSourceNodeId] |= (1LL << eventChannel);
        }
        else
        {
            ttlState[eventSourceNodeId] &= ~(1LL << eventChannel);
        }

        if (eventSourceNodeId == subprocessorToDraw)
        {
            const int chan          = numChannelsInSubprocessor[eventSourceNodeId];
            const int index         = (displayBufferIndex[chan] + eventTime) % displayBuffer->getNumSamples();
            const int samplesLeft   = displayBuffer->getNumSamples() - index;
            const int nSamples      = getNumSourceSamples(eventSourceNodeId) - eventTime;

            if (nSamples < samplesLeft)
            {
                displayBuffer->copyFrom(chan,                                 // destChannel
                                        index,                                // destStartSample
                                        arrayOfOnes,                          // source
                                        nSamples,                             // numSamples
                                        float(ttlState[eventSourceNodeId]));  // gain
            }
            else
            {
                int extraSamples = nSamples - samplesLeft;

                displayBuffer->copyFrom(chan,                                 // destChannel
                                        index,                                // destStartSample
                                        arrayOfOnes,                          // source
                                        samplesLeft,                          // numSamples
                                        float(ttlState[eventSourceNodeId]));  // gain

                displayBuffer->copyFrom(chan,                                 // destChannel
                                        0,                                    // destStartSample
                                        arrayOfOnes,                          // source
                                        extraSamples,                         // numSamples
                                        float(ttlState[eventSourceNodeId]));  // gain
            }
        }

        //         std::cout << "Received event from " << eventSourceNodeId
        //                   << " on channel " << eventChannel
        //                   << " with value " << eventId
        //                   << " at timestamp " << event.getTimeStamp() << std::endl;
    }
}

void LfpDisplayNode::handleSpike(const SpikeChannel* spikeInfo, const MidiMessage& event, int samplePosition) {
    SpikeEventPtr newSpike = SpikeEvent::deserializeFromMessage(event, spikeInfo); // a smart pointer
    if (!newSpike) return;
    
    int electrodeNum = getSpikeChannelIndex(newSpike);
    int timeStamp = newSpike->getTimestamp();

    Electrode* e = electrodes[electrodeNum];

    e->mostRecentSpikes.set(e->currentSpikeIndex, newSpike.release());
    e->currentSpikeIndex++;

    std::cout << "Spike saved: " << e->mostRecentSpikes.size() << std::endl;


    //std::cout << "Spike received: (" << electrodeNum << ":" <<timeStamp << ") "<< std::endl;


    
}


void LfpDisplayNode::initializeEventChannels()
{

	//std::cout << "Initializing events..." << std::endl;

    const int chan          = numChannelsInSubprocessor[subprocessorToDraw];
    const int index         = displayBufferIndex[chan];
    const int samplesLeft   = displayBuffer->getNumSamples() - index;
	const int nSamples      = getNumSourceSamples(subprocessorToDraw);

	//std::cout << chan << " " << index << " " << samplesLeft << " " << nSamples << std::endl;
        
    if (nSamples < samplesLeft)
    {

        displayBuffer->copyFrom (chan,                                      // destChannel
                                 index,                                     // destStartSample
                                 arrayOfOnes,                               // source
                                 nSamples,                                  // numSamples
                                 float (ttlState[subprocessorToDraw]));     // gain
    }
    else
    {
        // At the end of buffer, copy the last bit, and rewrite from the beginning
        int extraSamples = nSamples - samplesLeft;

        displayBuffer->copyFrom (chan,                                      // destChannel
                                 index,                                     // destStartSample
                                 arrayOfOnes,                               // source
                                 samplesLeft,                               // numSamples
                                 float (ttlState[subprocessorToDraw]));     // gain

        displayBuffer->copyFrom (chan,                                      // destChannel
                                 0,                                         // destStartSample
                                 arrayOfOnes,                               // source
                                 extraSamples,                              // numSamples
                                 float (ttlState[subprocessorToDraw]));     // gain
    }
}

void LfpDisplayNode::finalizeEventChannels()
{
    const int chan          = numChannelsInSubprocessor[subprocessorToDraw];
    const int index         = displayBufferIndex[chan];
    const int samplesLeft   = displayBuffer->getNumSamples() - index;
    const int nSamples      = getNumSourceSamples(subprocessorToDraw);
        
    int newIdx = 0;
        
    if (nSamples < samplesLeft)
    {
        newIdx = index + nSamples;
    }
    else
    {
        newIdx = nSamples - samplesLeft;
    }
        
    displayBufferIndex.set(chan, newIdx);
}


void LfpDisplayNode::process (AudioSampleBuffer& buffer)
{
    // 1. place any new samples into the displayBuffer
    //std::cout << "Display node sample count: " << nSamples << std::endl; ///buffer.getNumSamples() << std::endl;

	if (true)
	{
		ScopedLock displayLock(displayMutex);

		if (true)
		{
			initializeEventChannels();
			checkForEvents(true); // see if we got any TTL events
			finalizeEventChannels();
		}

		if (true)
		{
			int channelIndex = -1;

			for (int chan = 0; chan < buffer.getNumChannels(); ++chan)
			{
				if (getDataSubprocId(chan) == subprocessorToDraw)
				{
					channelIndex++;
					const int samplesLeft = displayBuffer->getNumSamples() - displayBufferIndex[channelIndex];
					const int nSamples = getNumSamples(chan);

					if (nSamples < samplesLeft)
					{
                        //Still some space in the display buffer, copy all
						displayBuffer->copyFrom(channelIndex,                      // destChannel
							displayBufferIndex[channelIndex],  // destStartSample
							buffer,                    // source
							chan,                      // source channel
							0,                         // source start sample
							nSamples);                 // numSamples

						displayBufferIndex.set(channelIndex, displayBufferIndex[channelIndex] + nSamples);
					}
					else
					{
                        // Copy the remaining, rewrite from the beginning

						const int extraSamples = nSamples - samplesLeft;

						displayBuffer->copyFrom(channelIndex,                      // destChannel
							displayBufferIndex[channelIndex],  // destStartSample
							buffer,                    // source
							chan,                      // source channel
							0,                         // source start sample
							samplesLeft);              // numSamples

						displayBuffer->copyFrom(channelIndex,                      // destChannel
							0,                         // destStartSample
							buffer,                    // source
							chan,                      // source channel
							samplesLeft,               // source start sample
							extraSamples);             // numSamples

						displayBufferIndex.set(channelIndex, extraSamples);
					}
				}
			}
		}
	}
}

