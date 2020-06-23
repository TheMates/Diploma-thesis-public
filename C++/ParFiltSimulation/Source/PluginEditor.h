/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"

//==============================================================================
/**
*/
class ParFiltSimulationAudioProcessorEditor : public AudioProcessorEditor, public Slider::Listener, public Button::Listener
{
public:
	ParFiltSimulationAudioProcessorEditor(ParFiltSimulationAudioProcessor&);
	~ParFiltSimulationAudioProcessorEditor();

	//==============================================================================
	void paint(Graphics&) override;
	void resized() override;

	// LISTENERS
	void sliderValueChanged(Slider* slider) override;
	void buttonClicked(Button* button) override;

private:
	// This reference is provided as a quick way for your editor to
	// access the processor object that created it.
	ParFiltSimulationAudioProcessor& processor;
	Slider gainSlider;
	
	std::unique_ptr<AudioProcessorValueTreeState::SliderAttachment> sliderAttach;

	std::unique_ptr<TextButton> ChooseBmBtn, ChooseAmBtn, ChooseFIRBtn, SetUpFilterBtn, BypassBtn, InfoBtn;
	std::unique_ptr<TextEditor> AmPath, BmPath, FIRPath;
	std::unique_ptr<Label> AmNFiltLab, BmNFiltLab, NFIRLab;
	std::unique_ptr<FileChooser> fileChooser;

	

	bool loadCoeffsFromCsv(const std::string &filePath,int ID);

	std::string cesta;

	int btnPressedID;

	void loadCoeffs(int ID);

	void checkReady();
	bool readyToSetupFilter = false;

	int nA = -1, nB = -1, nFIR = 0;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ParFiltSimulationAudioProcessorEditor)
};
