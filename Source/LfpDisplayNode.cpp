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
    spikeBuffer = new AudioSampleBuffer(8, 100);

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

void LfpViewer::LfpDisplayNode::updateSpikeElectrodeInfo()
{
    // Find out the mapping between electrode and channels
    uint16 electrodeNo = 0;
    for (int i = 0; i < this->electrodes.size(); i++) {
        // mark down which electrode correpond to which channel
        Array<uint16> electrodeMapping;
        for (int j = 0; j < electrodes[i]->numChannels; j++) {
            electrodeMapping.add(electrodeNo);

            printf("Channel %d belongs to electrode %d\n", electrodeNo, i);
            electrodeNo += 1;

        }
        electrode2channel.add(electrodeMapping);
    }
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

int64 LfpViewer::LfpDisplayNode::getDisplayBufferStartTimestamp()
{
    //return the timestasmp of the beginning of the display buffer
    return displayBufferStartTimestamp;
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
    totalSamples = nSamples;
	int nInputs = getNumSubprocessorChannels();

	std::cout << "Resizing buffer. Samples: " << nSamples << ", Inputs: " << nInputs << std::endl;

	if (nSamples > 0 && nInputs > 0)
	{
        //Note: this buffer length is different from the timebase used to show the signal in LfpDisplayCanvas
		abstractFifo.setTotalSize(nSamples);
		displayBuffer->setSize(nInputs + 1, nSamples); // add extra channel for TTLs
		displayBuffer->clear();

		displayBufferIndex.clear();
		displayBufferIndex.insertMultiple(0, 0, nInputs + 1);

        displayBufferStartTimestamp = 0;
   /*     displayBufferStartTimestamp.clear();
        displayBufferStartTimestamp.insertMultiple(0, 0, nInputs + 1);*/

        //initilaize the spike buffer
        
        this->spikeBuffer->setSize(nInputs+1, nSamples);
        this->spikeBuffer.get()->clear();
        printf("LfpDisplaynode electrode size set to %d \n", this->spikeBuffer.get()->getNumChannels());
        

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

    /*
    We have a spikeBuffer that has the same size as the displayBuffer. Upon receiving a spike, we mark in the spikeBuffer the location
    of the spike. When displayCanvas display the signal, it will copy displayBuffer into screenBuffer,
    it also has another screenSpikeBuffer that downsample the spikeBuffer into the same size as screenBuffer.
    The element inside screenSpikeBUffer indicate whether there is a spike in that pixel. This information can be used
    by channelDisplay to draw the spike raster line

    */

    // When the pointer go out of scope, it will be destoryed
    SpikeEventPtr newSpike = SpikeEvent::deserializeFromMessage(event, spikeInfo); // a smart pointer
    if (!newSpike) return;

    auto timeStamp = newSpike->getTimestamp();

    this->spikesLeftOver.add(newSpike.release());

    // mark in the spike buffer
    auto relativeTime = timeStamp - this->displayBufferStartTimestamp;
    if (relativeTime < totalSamples) {
        // Since the latest spikes has timesatmp within the display buffer range, then all previous
        // spikes should also be within range
        // pop left over spikes 

        for (SpikeEventPtr& spike : spikesLeftOver) {
            int electrodeNum = getSpikeChannelIndex(spike);
            auto timeStamp = spike->getTimestamp();

            auto relativeTime = timeStamp - this->displayBufferStartTimestamp;

            jassert(relativeTime < totalSamples);

            for (uint16& chanNo : this->electrode2channel[electrodeNum]) {
                if (relativeTime>0)
                    this->spikeBuffer->setSample(chanNo, relativeTime, 1);
                    //printf("Marking spikes at electrode %d time %d\n", electrodeNum, relativeTime);
            }
 
        }

        spikesLeftOver.clear();
        
    }

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
    // 2. if reached the end of the displaybuffer, then wrap around to the beginning
    // 3. displayBufferIndex contains the last write index location in displayBuffer
    // 4. displayBufferStartTimestamp contains the timestamp at the beginning of the buffer
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

        if (this->electrode2channel.size() == 0)
            updateSpikeElectrodeInfo();

       
		if (true)
		{
			int channelIndex = -1;

			for (int chan = 0; chan < buffer.getNumChannels(); ++chan)
			{
				if (getDataSubprocId(chan) == subprocessorToDraw)
				{
					channelIndex++;
					const int samplesLeft = displayBuffer->getNumSamples() - displayBufferIndex[channelIndex];
					const int nSamples = getNumSamples(chan); // number of sample in the input buffer

					if (nSamples < samplesLeft) // sampleLeft is the samples left in the display buffer
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

                        //Record the timestamp at the beginning of buffer
                        //displayBufferStartTimestamp.set(chan,getTimestamp(0) + samplesLeft);
                        displayBufferStartTimestamp = getTimestamp(0) + samplesLeft;
                        this->spikeBuffer->clear(); //clear the spike buffer
                        if (chan == 0)
                            std::cout << "displayBufferStartTimestamp: " << displayBufferStartTimestamp << std::endl;
                        
					}
				}
			}
		}
	}
}

