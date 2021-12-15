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

#ifndef __LFPDISPLAYNODE_H_Alpha__
#define __LFPDISPLAYNODE_H_Alpha__

#include <ProcessorHeaders.h>
#include "LfpDisplayEditor.h"

#include <map>

class DataViewport;

namespace LfpViewer
{

	/**

	  Holds data in a displayBuffer to be used by the LfpDisplayCanvas
	  for rendering continuous data streams.

	  @see GenericProcessor, LfpDisplayEditor, LfpDisplayCanvas

	*/

	struct Electrode
	{
		String name;

		int numChannels;
		int recordIndex;
		int currentSpikeIndex;

		Array<float> displayThresholds;
		Array<float> detectorThresholds;

		OwnedArray<SpikeEvent> mostRecentSpikes;

		float bitVolts;

	};


	struct SimpleElectrode
	{
		String name;

		int numChannels;
		int prePeakSamples, postPeakSamples;
		int lastBufferIndex;
		bool isMonitored;
		int electrodeID;
		int sourceNodeId;
		int recordIndex;
		int currentSpikeIndex;


		Array<float> displayThresholds;
		Array<float> detectorThresholds;

		HeapBlock<int> channels;
		HeapBlock<double> thresholds;
		HeapBlock<bool> isActive;

		float bitVolts;

	};

	class LfpDisplayNode : public GenericProcessor
	{
	public:
		LfpDisplayNode();
		~LfpDisplayNode();

		AudioProcessorEditor* createEditor() override;

		void process(AudioSampleBuffer& buffer) override;

		void setParameter(int parameterIndex, float newValue) override;

		void updateSettings() override;

		bool enable()   override;
		bool disable()  override;

		void handleEvent(const EventChannel* eventInfo, const MidiMessage& event, int samplePosition = 0) override;
		void handleSpike(const SpikeChannel* spikeInfo, const MidiMessage& event, int samplePosition) override;


		AudioSampleBuffer* getDisplayBufferAddress() const { return displayBuffer; }
		AudioSampleBuffer* getSpikeBufferAddress() const { return spikeBuffer; }

		int getDisplayBufferIndex(int chan) const { return displayBufferIndex[chan]; }

		CriticalSection* getMutex() { return &displayMutex; }

		void setSubprocessor(uint32 sp);
		uint32 getSubprocessor() const;

		int getNumSubprocessorChannels();

		float getSubprocessorSampleRate(uint32 subprocId);

		uint32 getDataSubprocId(int chan) const;

		OwnedArray<SimpleElectrode>* getElectrodes();

		int64 getDisplayBufferStartTimestamp();

		StringArray electrodeTypes;
		StringArray detectionMethod; //method of detection
		StringArray detectionSign; //direction of detection
		String curDetectionMethod;
		String curDetectionSign;

		//Spike detection related
		void computeMedianThreshold(SimpleElectrode* electrode, int nSamples, int& sample_counter, std::vector<std::vector<float>>& dyn_thresholds, int& window_number);
		void computeSimpleThreshold(SimpleElectrode* electrode, std::vector<std::vector<float>>& thresholds);
		int getNumChannels(int index);
		double getChannelThreshold(int electrodeNum, int channelNum);
		void setChannelActive(int electrodeIndex, int subChannel, bool active);
		void setEnableDetection(bool isEnable);
		void setChannelThreshold(int electrodeNum, int channelNum, float thresh);
		bool removeElectrode(int index);
		bool addElectrode(int nChans, int electrodeID);
		void resetElectrode(SimpleElectrode* e);
		bool isChannelActive(int electrodeIndex, int channelNum);
		StringArray getElectrodeNames();
		float getDefaultThreshold();
		int getChannel(int index, int i);
		SimpleElectrode* getActiveElectrode();


	private:
		void initializeEventChannels();
		void finalizeEventChannels();

		OwnedArray<SimpleElectrode> electrodes;

		ScopedPointer<AudioSampleBuffer> displayBuffer;

		Array<int> displayBufferIndex;
		Array<uint32> eventSourceNodes;

		float displayGain; //
		float bufferLength; // 
		int64 displayBufferStartTimestamp; //timestamp of the beginning of display buffer, in no. of samples, not actual time

		AbstractFifo abstractFifo;

		int64 bufferTimestamp;
		std::map<uint32, uint64> ttlState;
		float* arrayOfOnes;
		int totalSamples;

		bool resizeBuffer();

		int numSubprocessors;
		uint32 subprocessorToDraw;
		std::map<uint32, int> numChannelsInSubprocessor;
		std::map<uint32, float> subprocessorSampleRate;

		CriticalSection displayMutex;

		static uint32 getEventSourceId(const EventChannel* event);
		static uint32 getChannelSourceId(const InfoObjectCommon* chan);


		Array<Array<uint16>> electrode2channel;
		void updateSpikeElectrodeInfo();

		ScopedPointer<AudioSampleBuffer> spikeBuffer; //whether a certain pixel in a certain pixel contains spikes
		Array<SpikeEventPtr> spikesLeftOver; //leftover spikes

		//////////////////////
		// Spike detection related functions
		Array<int> electrodeCounter;
		bool isEnableDetection = false;
		int uniqueID;
		int currentElectrode;
		int currentChannelIndex;
		int currentIndex;


		JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LfpDisplayNode);

	};
};



#endif  // __LFPDISPLAYNODE_H_Alpha__
