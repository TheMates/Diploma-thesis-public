/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
ParFiltSimulationAudioProcessor::ParFiltSimulationAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", AudioChannelSet::stereo(), true)
                     #endif
                       ),
	gainValue(0.),
	parameters(*this, nullptr )
#endif
{

	filters.emplace_back(std::make_unique<FilterBank>(Am, Am, Bm, Bm, 0, FIR, 0, 0.f));
	filters.emplace_back(std::make_unique<FilterBank>(Am, Am, Bm, Bm, 0, FIR, 0, 0.f));

	NormalisableRange<float> gainRange(-48.f,0.f);
	parameters.createAndAddParameter("gain", "Gain", "Gain", gainRange, 0.f, nullptr, nullptr);
	parameters.state = ValueTree("savedParams");

}

ParFiltSimulationAudioProcessor::~ParFiltSimulationAudioProcessor()
{
}

//==============================================================================
const String ParFiltSimulationAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool ParFiltSimulationAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool ParFiltSimulationAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool ParFiltSimulationAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double ParFiltSimulationAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int ParFiltSimulationAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int ParFiltSimulationAudioProcessor::getCurrentProgram()
{
    return 0;
}

void ParFiltSimulationAudioProcessor::setCurrentProgram (int index)
{
}

const String ParFiltSimulationAudioProcessor::getProgramName (int index)
{
    return {};
}

void ParFiltSimulationAudioProcessor::changeProgramName (int index, const String& newName)
{
}

//==============================================================================
void ParFiltSimulationAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    // Use this method as the place to do any pre-playback
    // initialisation that you need..
}

void ParFiltSimulationAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool ParFiltSimulationAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    if (layouts.getMainOutputChannelSet() != AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

void ParFiltSimulationAudioProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
    ScopedNoDenormals noDenormals;
    auto nInputChan  = getTotalNumInputChannels();
    auto nOutputChan = getTotalNumOutputChannels();

    //clears the output if there is more inputs than outputs
    for (auto i = nInputChan; i < nOutputChan; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());

    for (int channel = 0; channel < nInputChan; ++channel)
    {
        auto* channelData = buffer.getWritePointer (channel);
		

		if (filterSet && !isBypassed)
		{
			for (auto sample = 0; sample < buffer.getNumSamples(); ++sample)
			{
				float processed = filters[channel]->processSample(buffer.getSample(channel, sample));
				channelData[sample] = processed * pow(10., gainValue/20.);
			}
		}
		else
		{
			for (auto sample = 0; sample < buffer.getNumSamples(); ++sample)
			{
				float s = buffer.getSample(channel, sample);
				channelData[sample] = buffer.getSample(channel, sample)* pow(10., gainValue / 20.);
			}
		}
    }
}

//==============================================================================
bool ParFiltSimulationAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* ParFiltSimulationAudioProcessor::createEditor()
{
	return new ParFiltSimulationAudioProcessorEditor (*this);
}

//==============================================================================
void ParFiltSimulationAudioProcessor::getStateInformation (MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.

	std::unique_ptr<XmlElement> xml( parameters.state.createXml() );
	copyXmlToBinary(*xml, destData);
}

void ParFiltSimulationAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
	std::unique_ptr<XmlElement> savedParams(getXmlFromBinary(data, sizeInBytes));
	if (savedParams != nullptr)
		if (savedParams->hasTagName(parameters.state.getType()))
			parameters.state = ValueTree::fromXml(*savedParams);

}

void ParFiltSimulationAudioProcessor::filterChanged()
{
	filterSet = false;
	filters[0]->Reset();
	filters[1]->Reset();
}

void ParFiltSimulationAudioProcessor::setUpFilter()
{
	filters[0].reset();
	filters[0] = std::make_unique<FilterBank>(Am, Am + nFilters, Bm, Bm + nFilters, nFilters, FIR, FIRsize, 1.0);
	filters[1].reset();
	filters[1] = std::make_unique<FilterBank>(Am, Am + nFilters, Bm, Bm + nFilters, nFilters, FIR, FIRsize, 1.0);
	filterSet = true;
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new ParFiltSimulationAudioProcessor();
}
