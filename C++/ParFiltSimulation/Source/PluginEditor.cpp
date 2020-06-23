/*
  ==============================================================================

	This file was auto-generated!

	It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"
#include "JuceUTF16toStdString.h"
#include "rapidcsv.h"
#include "Windows.h"

using namespace std;

//==============================================================================
ParFiltSimulationAudioProcessorEditor::ParFiltSimulationAudioProcessorEditor(ParFiltSimulationAudioProcessor& p)
	: AudioProcessorEditor(&p), processor(p)
{
	// Make sure that before the constructor has finished, you've set the
	// editor's size to whatever you need it to be.
	

	gainSlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
	gainSlider.setTextBoxStyle(Slider::TextBoxBelow, false, 75, 25);
	gainSlider.setRange(-48., 0.);
	gainSlider.addListener(this);
	//gainSlider.setValue(-6);
	addAndMakeVisible(gainSlider);

	sliderAttach = std::make_unique<AudioProcessorValueTreeState::SliderAttachment>(processor.parameters,"gain",gainSlider);


	ChooseAmBtn = std::make_unique<TextButton>();
	ChooseAmBtn->setButtonText("Load Am coeffs");
	ChooseAmBtn->addListener(this);
	ChooseAmBtn->setComponentID("0");
	addAndMakeVisible(ChooseAmBtn.get());

	ChooseBmBtn = std::make_unique<TextButton>();
	ChooseBmBtn->setButtonText("Load Bm coeffs");
	ChooseBmBtn->addListener(this);
	ChooseBmBtn->setComponentID("1");
	addAndMakeVisible(ChooseBmBtn.get());

	ChooseFIRBtn = std::make_unique<TextButton>();
	ChooseFIRBtn->setButtonText("Load FIR coeffs");
	ChooseFIRBtn->addListener(this);
	ChooseFIRBtn->setComponentID("2");
	addAndMakeVisible(ChooseFIRBtn.get());

	SetUpFilterBtn = std::make_unique<TextButton>();
	SetUpFilterBtn->setButtonText("Set Filter");
	SetUpFilterBtn->setColour(TextButton::ColourIds::buttonColourId, Colours::red);
	SetUpFilterBtn->addListener(this);
	addAndMakeVisible(SetUpFilterBtn.get());

	BypassBtn = std::make_unique<TextButton>();
	BypassBtn->setButtonText("Bypass");
	BypassBtn->setColour(TextButton::ColourIds::buttonColourId, Colours::grey);
	BypassBtn->addListener(this);
	BypassBtn->setComponentID("8");
	addAndMakeVisible(BypassBtn.get());

	InfoBtn = std::make_unique<TextButton>();
	InfoBtn->setButtonText("i");
	InfoBtn->setColour(TextButton::ColourIds::buttonColourId, Colours::lightsteelblue);
	InfoBtn->addListener(this);
	InfoBtn->setComponentID("9");
	addAndMakeVisible(InfoBtn.get());

	AmPath = std::make_unique<TextEditor>();
	AmPath->setEnabled(false);
	addAndMakeVisible(AmPath.get());

	BmPath = std::make_unique<TextEditor>();
	BmPath->setEnabled(false);
	addAndMakeVisible(BmPath.get());

	FIRPath = std::make_unique<TextEditor>();
	FIRPath->setEnabled(false);
	addAndMakeVisible(FIRPath.get());

	AmNFiltLab = std::make_unique<Label>();
	addAndMakeVisible(AmNFiltLab.get());

	BmNFiltLab = std::make_unique<Label>();
	addAndMakeVisible(BmNFiltLab.get());
	
	NFIRLab = std::make_unique<Label>();
	addAndMakeVisible(NFIRLab.get());
	
	setSize(400, 160);



}

ParFiltSimulationAudioProcessorEditor::~ParFiltSimulationAudioProcessorEditor()
{
}

//==============================================================================
void ParFiltSimulationAudioProcessorEditor::paint(Graphics& g)
{
	// (Our component is opaque, so we must completely fill the background with a solid colour)
	g.fillAll(getLookAndFeel().findColour(ResizableWindow::backgroundColourId));


	//g.setColour (Colours::white);
	//g.setFont (15.0f);
	//g.drawFittedText ("Hello World!", getLocalBounds(), Justification::centred, 1);
}

void ParFiltSimulationAudioProcessorEditor::resized()
{
	// This is generally where you'll want to lay out the positions of any
	// subcomponents in your editor..

	auto area = getLocalBounds();
	ChooseAmBtn->setBounds(10, 10, 120, 25);
	ChooseBmBtn->setBounds(10, 37, 120, 25);
	ChooseFIRBtn->setBounds(10, 64, 120, 25);
	SetUpFilterBtn->setBounds(10, 95, 120, 25);
	BypassBtn->setBounds(40, 125, 90, 25);
	InfoBtn->setBounds(10, 125, 25, 25);

	AmPath->setBounds(160, 10, 150, 25);
	BmPath->setBounds(160, 37, 150, 25);
	FIRPath->setBounds(160, 64, 150, 25);

	AmNFiltLab->setBounds(130, 10, 30, 25);
	BmNFiltLab->setBounds(130, 37, 30, 25);
	NFIRLab->setBounds(130, 64, 20, 25);

	gainSlider.setBounds(area.getWidth() - 75, 0, 75, getHeight());
}

void ParFiltSimulationAudioProcessorEditor::sliderValueChanged(Slider* slider)
{
	if (slider == &gainSlider)
	{
		processor.gainValue = slider->getValue();
	}
}

void ParFiltSimulationAudioProcessorEditor::buttonClicked(Button* button)
{
	if (button == SetUpFilterBtn.get())
	{
		if (readyToSetupFilter)
		{
			processor.setUpFilter();
			SetUpFilterBtn->setColour(TextButton::ColourIds::buttonColourId, Colours::green);
		}

		return;
	}
	if(button == BypassBtn.get())
	{
		processor.isBypassed = !processor.isBypassed;
		if(processor.isBypassed)
			BypassBtn->setColour(TextButton::ColourIds::buttonColourId, Colours::red);
		else
			BypassBtn->setColour(TextButton::ColourIds::buttonColourId, Colours::grey);
		return;
	}
	if(button == InfoBtn.get())
	{
		AlertWindow::showMessageBoxAsync(AlertWindow::InfoIcon,
			"About this plug-in",
			CharPointer_UTF8("Plug-in simulates frequency response by parallel filter of second-order sections.\n"
				"\n"
				"Load csv documents comtaining coefficients of these sections. \n"
				"Am - denominator coefficients, 3 column vector.\n"
				"Bm - numerator coefficients, 2 column vector.\n"
				"FIR - FIR coefficients, not necessary, 1 column vector, or scalar.\n"
				"\n"
				"For information about parallel filters see:\n"
				"http://home.mit.bme.hu/~bank/parfilt/\n"
				"\n"
				"Plug-in by Matou\xc5\xa1 Vrb\xc3\xadk\n"
				"matousvrbik@gmail.com")
		);
	}

	btnPressedID = std::stoi(button->getComponentID().toStdString());

	switch (btnPressedID)
	{
	case 0:	//Load Am		
	case 1:	//Load Bm			
	case 2:	//Load FIR
		fileChooser.reset(new FileChooser("Load csv file", File::getCurrentWorkingDirectory(),
			"*.csv", true));
		fileChooser->launchAsync(FileBrowserComponent::openMode | FileBrowserComponent::canSelectFiles,
			[&](const FileChooser& chooser)
		{
			String chosen;
			auto results = chooser.getURLResults();

			for (auto result : results)
				chosen << (result.isLocalFile() ? result.getLocalFile().getFullPathName()
					: result.toString(false));

			cesta = JuceUTF16toStdString::Converter::Convert(chosen.toUTF16());
			//AlertWindow::showMessageBoxAsync(AlertWindow::InfoIcon,
			//	"File Chooser...",
			//	"You picked: " + chosen);
			loadCoeffs(btnPressedID);
			processor.filterChanged();
		});
		break;

	default: break;
	}
}

bool ParFiltSimulationAudioProcessorEditor::loadCoeffsFromCsv(const std::string &filePath, int ID)
{
	float* coeffs = nullptr;
	Label* lab;
	TextEditor* textPath;
	int* nFilt;
	switch (ID)
	{
	case 0: lab = AmNFiltLab.get(); textPath = AmPath.get(); nFilt = &nA; break;
	case 1: lab = BmNFiltLab.get(); textPath = BmPath.get(); nFilt = &nB; break;
	case 2: lab = NFIRLab.get(); textPath = FIRPath.get(); nFilt = &nFIR; break;
	default: break;
	}
	try
	{
		if (filePath == "")
			throw exception("0");
		//choose region dependent separator
		char sep[2];
		GetLocaleInfoA(LOCALE_USER_DEFAULT, LOCALE_SLIST, sep, sizeof(sep));
		rapidcsv::Document doc(filePath, rapidcsv::LabelParams(-1, -1),rapidcsv::SeparatorParams(sep[0]));			//so first column and row are counted as data, not headers
		const int rows = doc.GetRowCount();
		const int cols = doc.GetColumnCount();

		int colOffset = cols == 3 ? 1 : 0;

		coeffs = new float[(cols - colOffset) * rows];	//allocate memory

		vector<float> temp;
		for (auto col = 0; col < cols - colOffset; ++col)
		{
			temp = doc.GetColumn<float>(col+ colOffset);
			std::copy(temp.begin(), temp.end(), coeffs + col*rows);
		}

		lab->setText(std::to_string(rows),NotificationType::dontSendNotification);
		*nFilt = rows;
		textPath->setText(filePath, false);
		switch (ID)
		{
		case 0: processor.Am = coeffs;  break;
		case 1: processor.Bm = coeffs;  break;
		case 2: processor.FIR = coeffs; break;
		default: break;
		}

	}
	catch (exception& e)
	{
		if (coeffs != nullptr)
			delete[] coeffs;
		coeffs = nullptr;
		lab->setText("", NotificationType::dontSendNotification);
		*nFilt = -1;
		textPath->setText("", false);
		
		switch (ID)
		{
		case 0: processor.Am = nullptr;  break;
		case 1: processor.Bm = nullptr;  break;
		case 2: processor.FIR = nullptr; *nFilt = 0; break;
		default: break;
		}
		if (strcmp(e.what(),"0") != 0)	//no file selected -> no message necessary
		{
			AlertWindow::showMessageBoxAsync(AlertWindow::WarningIcon,
				"Loading error!",
				e.what());
		}
		return false;
	}
	return true;
}

void ParFiltSimulationAudioProcessorEditor::loadCoeffs(int ID)
{
	switch (ID)
	{
	case 0:	//Am
		if(processor.Am != nullptr) 
			delete[] processor.Am;
		processor.Am = nullptr;
		processor.AmSet = loadCoeffsFromCsv(cesta, ID);
		break;
	case 1:	//Bm
		if (processor.Bm != nullptr)
			delete[] processor.Bm;
		processor.Bm = nullptr;
		processor.BmSet = loadCoeffsFromCsv(cesta, ID);
		break;

	case 2:	//FIR
		if (processor.FIR != nullptr)
			delete[] processor.FIR;
		processor.FIR = nullptr;
		processor.FIRSet = loadCoeffsFromCsv(cesta, ID);
		break;

	default:break;
	}
	checkReady();
}

void ParFiltSimulationAudioProcessorEditor::checkReady()
{
	if(nA>0 && nA == nB && nFIR>=0)	//ready to set up
	{
		SetUpFilterBtn->setColour(TextButton::ColourIds::buttonColourId, Colours::blue);
		readyToSetupFilter = true;
		processor.nFilters = nA;
		processor.FIRsize = nFIR;
	}
	else
	{
		SetUpFilterBtn->setColour(TextButton::ColourIds::buttonColourId, Colours::red);
		readyToSetupFilter = false;
		processor.filterChanged();
	}
}



