/*
    ------------------------------------------------------------------

    This file is part of the Open Ephys GUI
    Copyright (C) 2013 Open Ephys

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

#ifndef __LFPDISPLAYEDITOR_H_Alpha__
#define __LFPDISPLAYEDITOR_H_Alpha__

#include <VisualizerEditorHeaders.h>
#include "LfpDisplayNode.h"
#include "LfpDisplayCanvas.h"

class Visualizer;

namespace LfpViewer {
    
class LfpDisplayNode;
class LfpDisplayCanvas;

/**

  User interface for the LfpDisplayNode sink.

  @see LfpDisplayNode, LfpDisplayCanvas

*/

class LfpDisplayEditor : public VisualizerEditor,
    public ComboBox::Listener,
    public Label::Listener
{
public:
    LfpDisplayEditor(GenericProcessor*, bool useDefaultParameterEditors);
    ~LfpDisplayEditor();

    /** Override the default VisualizerEditor behavior slightly, only for
        initialization 
     */
    void buttonClicked(Button* button) override;
    // not really being used (yet) ...
    void buttonEvent(Button* button);
    /** Respond to user's subprocessor sample rate selection */
    void comboBoxChanged(ComboBox *cb);

    /** Called by the base class VisualizerEditor to display the canvas
        when the user chooses to display one
     
        @see VisualizerEditor::buttonClicked
     */
    Visualizer* createNewCanvas();

	void startAcquisition();
	void stopAcquisition();

	void saveVisualizerParameters(XmlElement* xml);
	void loadVisualizerParameters(XmlElement* xml);
    
    /** Handle the state and options within the subprocessor sample rate
        selection combobox 
     */
    void updateSubprocessorSelectorOptions();

private:
    
    SortedSet<uint32> inputSubprocessors;
    
    LfpDisplayNode* lfpProcessor;

    // label and combobox for subprocessor selection
    ScopedPointer<Label> subprocessorSelectionLabel;
    ScopedPointer<ComboBox> subprocessorSelection;
    
    ScopedPointer<Label> subprocessorSampleRateLabel;
    
    bool hasNoInputs;

	int defaultSubprocessor;

    // TODO: use managed ptr to avoid memory leak
    ComboBox* electrodeTypes;
    ComboBox* electrodeList;
    Label* numElectrodes;
    Label* thresholdLabel;
    TriangleButton* upButton;
    TriangleButton* downButton;
    UtilityButton* plusButton;
    ToggleButton* enableBtn;

    ThresholdSlider* thresholdSlider;

    OwnedArray<ElectrodeButton> electrodeButtons;
    Array<ElectrodeEditorButton*> electrodeEditorButtons;

    ComboBox* detectionMethod;
    ComboBox* detectionSign;
    Label* detectionMethodLabel;
    Label* thresholdTextBox;
    TextButton* setAllThresholdBtn; //button to set all threshold

    void editElectrode(int index, int chan, int newChan);

    int lastId;
    bool isPlural;

    void labelTextChanged(Label* label);
    Font font;

    bool addElectrode(int nChans, int electrodeID=0);
    void refreshElectrodeList();
    void drawElectrodeButtons(int ID);
    void removeElectrode(int index);

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LfpDisplayEditor);

};
};
#endif  // __LFPDISPLAYEDITOR_H_Alpha__
