package simulation.software.codebase;

import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.beans.property.DoubleProperty;
import javafx.beans.value.ChangeListener;
import javafx.fxml.FXML;
import javafx.scene.control.*;
import javafx.scene.layout.VBox;
import javafx.stage.FileChooser;
import javafx.util.Duration;
import java.io.File;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.logging.Logger;

public class FilterController {
    private static final Logger LOGGER = Logger.getLogger(FilterController.class.getName());

    @FXML private TabPane mainTabPane;
    @FXML private ComboBox<String> filterType;
    @FXML private ComboBox<String> designType;
    @FXML private ComboBox<String> designMethod;
    @FXML private TextField orderField;
    @FXML private CheckBox autoOrderCheck;
    @FXML private Slider orderSlider;
    @FXML private TextField passbandRippleField;
    @FXML private Slider passbandRippleSlider;
    @FXML private TextField stopbandAttenuationField;
    @FXML private Slider stopbandAttenuationSlider;
    @FXML private TextField cutoffFrequency1Field;
    @FXML private Slider cutoffFrequency1Slider;
    @FXML private TextField cutoffFrequency2Field;
    @FXML private Slider cutoffFrequency2Slider;
    @FXML private Label cutoff2Label;
    @FXML private TextField passbandEdgeField;
    @FXML private Slider passbandEdgeSlider;
    @FXML private TextField stopbandEdgeField;
    @FXML private Slider stopbandEdgeSlider;
    @FXML private ComboBox<String> testSignalType;
    @FXML private TextField testSignalFreqField;
    @FXML private Slider testSignalFreqSlider;
    @FXML private TextField testSignalAmpField;
    @FXML private Slider testSignalAmpSlider;
    @FXML private Label stabilityLabel;
    @FXML private VBox magDbContainer;
    @FXML private VBox magLinearContainer;
    @FXML private VBox phaseContainer;
    @FXML private VBox groupDelayContainer;
    @FXML private VBox poleZeroContainer;
    @FXML private VBox impulseContainer;
    @FXML private VBox stepContainer;
    @FXML private VBox signalContainer;
    @FXML private Button saveButton;
    @FXML private Button startSimButton;
    @FXML private Button applySignalButton;
    @FXML private Button resetButton;
    @FXML private TextArea designOutputArea;
    @FXML private Label statusLabel;

    private FilterCharts charts;
    private ChebyshevFilter currentFilter;
    private Timeline debounceTimeline;
    private Timeline sweepTimeline;
    private AtomicBoolean isSweeping = new AtomicBoolean(false);

    @FXML
    public void initialize() {
        filterType.getItems().addAll("Type I", "Type II");
        filterType.setValue("Type I");
        designType.getItems().addAll("Low-pass", "High-pass", "Band-pass", "Band-stop");
        designType.setValue("Low-pass");
        designMethod.getItems().addAll("Transfer Function", "ZPK", "SOS");
        designMethod.setValue("ZPK");
        testSignalType.getItems().addAll("Sine", "Square", "Triangle", "Sawtooth");
        testSignalType.setValue("Sine");
        charts = new FilterCharts();
        magDbContainer.getChildren().add(charts.getMagDbPane());
        magLinearContainer.getChildren().add(charts.getMagLinearPane());
        phaseContainer.getChildren().add(charts.getPhasePane());
        groupDelayContainer.getChildren().add(charts.getGroupDelayPane());
        poleZeroContainer.getChildren().add(charts.getPoleZeroPane());
        impulseContainer.getChildren().add(charts.getImpulsePane());
        stepContainer.getChildren().add(charts.getStepPane());
        signalContainer.getChildren().add(charts.getSignalPane());

        // Debounce updates
        debounceTimeline = new Timeline(new KeyFrame(Duration.millis(300), e -> updateFilter()));
        debounceTimeline.setCycleCount(1);

        // Bind sliders to text fields
        bindSliderToTextField(orderSlider.valueProperty(), orderField, "%.0f");
        bindSliderToTextField(passbandRippleSlider.valueProperty(), passbandRippleField, "%.2f");
        bindSliderToTextField(stopbandAttenuationSlider.valueProperty(), stopbandAttenuationField, "%.2f");
        bindSliderToTextField(cutoffFrequency1Slider.valueProperty(), cutoffFrequency1Field, "%.2f");
        bindSliderToTextField(cutoffFrequency2Slider.valueProperty(), cutoffFrequency2Field, "%.2f");
        bindSliderToTextField(passbandEdgeSlider.valueProperty(), passbandEdgeField, "%.2f");
        bindSliderToTextField(stopbandEdgeSlider.valueProperty(), stopbandEdgeField, "%.2f");
        bindSliderToTextField(testSignalFreqSlider.valueProperty(), testSignalFreqField, "%.2f");
        bindSliderToTextField(testSignalAmpSlider.valueProperty(), testSignalAmpField, "%.2f");

        // Initialize UI
        designType.valueProperty().addListener((obs, oldVal, newVal) -> {
            boolean isBand = newVal.equals("Band-pass") || newVal.equals("Band-stop");
            cutoff2Label.setVisible(isBand);
            cutoffFrequency2Field.setVisible(isBand);
            cutoffFrequency2Slider.setVisible(isBand);
            triggerUpdate();
        });

        autoOrderCheck.selectedProperty().addListener((obs, oldVal, newVal) -> {
            orderField.setDisable(newVal);
            orderSlider.setDisable(newVal);
            triggerUpdate();
        });

        // Slider listeners
        orderSlider.valueProperty().addListener((obs, oldVal, newVal) -> triggerUpdate());
        passbandRippleSlider.valueProperty().addListener((obs, oldVal, newVal) -> triggerUpdate());
        stopbandAttenuationSlider.valueProperty().addListener((obs, oldVal, newVal) -> triggerUpdate());
        cutoffFrequency1Slider.valueProperty().addListener((obs, oldVal, newVal) -> triggerUpdate());
        cutoffFrequency2Slider.valueProperty().addListener((obs, oldVal, newVal) -> triggerUpdate());
        passbandEdgeSlider.valueProperty().addListener((obs, oldVal, newVal) -> triggerUpdate());
        stopbandEdgeSlider.valueProperty().addListener((obs, oldVal, newVal) -> triggerUpdate());
        testSignalFreqSlider.valueProperty().addListener((obs, oldVal, newVal) -> triggerUpdate());
        testSignalAmpSlider.valueProperty().addListener((obs, oldVal, newVal) -> triggerUpdate());

        // Text field listeners
        ChangeListener<String> updateListener = (obs, oldVal, newVal) -> triggerUpdate();
        for (TextField field : new TextField[]{orderField, passbandRippleField, stopbandAttenuationField,
                cutoffFrequency1Field, cutoffFrequency2Field, passbandEdgeField, stopbandEdgeField,
                testSignalFreqField, testSignalAmpField}) {
            field.textProperty().addListener(updateListener);
        }
        filterType.valueProperty().addListener(updateListener);
        designMethod.valueProperty().addListener((obs, oldVal, newVal) -> updateDesignOutput());
        testSignalType.valueProperty().addListener(updateListener);

        // Button actions
        saveButton.setOnAction(e -> saveConfig());
        startSimButton.setOnAction(e -> toggleSweep());
        applySignalButton.setOnAction(e -> applyTestSignal());
        resetButton.setOnAction(e -> resetFields());

        // Tooltips
        filterType.setTooltip(new Tooltip("Select Chebyshev filter type"));
        designType.setTooltip(new Tooltip("Select filter design type"));
        orderSlider.setTooltip(new Tooltip("Adjust filter order"));
        cutoffFrequency1Slider.setTooltip(new Tooltip("Adjust first cutoff frequency"));
        cutoffFrequency2Slider.setTooltip(new Tooltip("Adjust second cutoff frequency"));
        passbandRippleSlider.setTooltip(new Tooltip("Adjust passband ripple"));
        stopbandAttenuationSlider.setTooltip(new Tooltip("Adjust stopband attenuation"));
        passbandEdgeSlider.setTooltip(new Tooltip("Adjust passband edge frequency"));
        stopbandEdgeSlider.setTooltip(new Tooltip("Adjust stopband edge frequency"));
        testSignalType.setTooltip(new Tooltip("Select test signal type"));
        testSignalFreqSlider.setTooltip(new Tooltip("Adjust test signal frequency"));
        testSignalAmpSlider.setTooltip(new Tooltip("Adjust test signal amplitude"));

        LOGGER.info("FilterController initialized");
        resetFields();
    }

    private void bindSliderToTextField(DoubleProperty sliderProperty, TextField textField, String format) {
        textField.textProperty().bindBidirectional(sliderProperty, new javafx.util.StringConverter<Number>() {
            @Override
            public String toString(Number value) {
                return String.format(format, value.doubleValue());
            }

            @Override
            public Number fromString(String string) {
                try {
                    return Double.parseDouble(string);
                } catch (NumberFormatException e) {
                    return sliderProperty.get();
                }
            }
        });
    }

    private void triggerUpdate() {
        statusLabel.setText("Computing...");
        debounceTimeline.stop();
        debounceTimeline.playFromStart();
        LOGGER.info("Triggered filter update");
    }

    private void updateFilter() {
        try {
            if (!validateInputs()) {
                showAlert("Invalid input. Please check all fields.");
                designOutputArea.setText("Invalid inputs");
                statusLabel.setText("Error: Invalid inputs");
                stabilityLabel.setText("Unknown");
                LOGGER.warning("Invalid inputs detected");
                return;
            }

            ChebyshevFilter.FilterType filterTypeVal = filterType.getValue().equals("Type I") ?
                    ChebyshevFilter.FilterType.TYPE_I : ChebyshevFilter.FilterType.TYPE_II;
            ChebyshevFilter.DesignType designTypeVal = ChebyshevFilter.DesignType.valueOf(
                    designType.getValue().replace("-", "_").toUpperCase());
            double passbandRipple = parseDouble(passbandRippleField.getText(), 1.0);
            double stopbandAttenuation = parseDouble(stopbandAttenuationField.getText(), 20.0);
            double cutoffFrequency1 = parseDouble(cutoffFrequency1Field.getText(), 1.0);
            double cutoffFrequency2 = parseDouble(cutoffFrequency2Field.getText(),
                    designTypeVal == ChebyshevFilter.DesignType.BAND_PASS ||
                            designTypeVal == ChebyshevFilter.DesignType.BAND_STOP ? 2.0 : cutoffFrequency1);
            double passbandEdge = parseDouble(passbandEdgeField.getText(), 1.0);
            double stopbandEdge = parseDouble(stopbandEdgeField.getText(), 2.0);

            int order;
            if (autoOrderCheck.isSelected()) {
                order = ChebyshevFilter.calculateMinOrder(passbandRipple, stopbandAttenuation, passbandEdge, stopbandEdge);
                orderField.setText(String.valueOf(order));
                orderSlider.setValue(order);
            } else {
                order = parseInt(orderField.getText(), 4);
            }

            currentFilter = new ChebyshevFilter(filterTypeVal, designTypeVal, order, passbandRipple,
                    stopbandAttenuation, cutoffFrequency1, cutoffFrequency2, passbandEdge, stopbandEdge);

            double startFreq = 0.1;
            double endFreq = designTypeVal == ChebyshevFilter.DesignType.BAND_PASS ||
                    designTypeVal == ChebyshevFilter.DesignType.BAND_STOP ? 10.0 : 5.0;
            double duration = 5.0;
            double dt = 0.01;
            double[] frequencies = currentFilter.getFrequencies(startFreq, endFreq, 1000);
            double[] magDb = currentFilter.getFrequencyResponse(startFreq, endFreq, 1000);
            double[] magLinear = currentFilter.getMagnitudeResponseLinear(startFreq, endFreq, 1000);
            double[] phase = currentFilter.getPhaseResponse(startFreq, endFreq, 1000);
            double[] groupDelay = currentFilter.getGroupDelay(startFreq, endFreq, 1000);
            double[] time = generateTimeVector(duration, dt);
            double[] impulse = currentFilter.getImpulseResponse(duration, dt);
            double[] step = currentFilter.getStepResponse(duration, dt);

            LOGGER.info("Updating plots for " + filterTypeVal + " " + designTypeVal);
            charts.updateFilterCharts(frequencies, magDb, magLinear, phase, groupDelay,
                    currentFilter.getPoles(), currentFilter.getZeros(), time, impulse, step,
                    filterTypeVal + " " + designTypeVal, true);
            updateDesignOutput();
            applyTestSignal();

            // Update stability
            boolean isStable = currentFilter.isStable();
            stabilityLabel.setText(isStable ? "Stable" : "Unstable");
            stabilityLabel.setStyle(isStable ? "-fx-text-fill: green;" : "-fx-text-fill: red;");
            statusLabel.setText("Plots updated, " + (isStable ? "stable" : "unstable") + " design");
            LOGGER.info("Stability check: " + (isStable ? "Stable" : "Unstable"));
        } catch (Exception e) {
            showAlert("Error: " + e.getMessage());
            designOutputArea.setText("Error computing filter");
            statusLabel.setText("Error: " + e.getMessage());
            stabilityLabel.setText("Unknown");
            LOGGER.severe("Error updating filter: " + e.getMessage());
        }
    }

    private void updateDesignOutput() {
        if (currentFilter == null) return;
        ChebyshevFilter.DesignMethod method = ChebyshevFilter.DesignMethod.valueOf(
                designMethod.getValue().replace(" ", "_").toUpperCase());
        designOutputArea.setText(currentFilter.getDesignRepresentation(method));
        LOGGER.info("Design output updated for method: " + method);
    }

    private boolean validateInputs() {
        try {
            for (TextField field : new TextField[]{passbandRippleField, stopbandAttenuationField, cutoffFrequency1Field, passbandEdgeField, stopbandEdgeField}) {
                if (!field.getText().isEmpty() && Double.parseDouble(field.getText()) <= 0) {
                    return false;
                }
            }
            if (cutoffFrequency2Field.isVisible() && !cutoffFrequency2Field.getText().isEmpty() && Double.parseDouble(cutoffFrequency2Field.getText()) <= 0) {
                return false;
            }
            if (!orderField.getText().isEmpty() && Integer.parseInt(orderField.getText()) <= 0) {
                return false;
            }
            if (cutoffFrequency2Field.isVisible() && !cutoffFrequency1Field.getText().isEmpty() && !cutoffFrequency2Field.getText().isEmpty()) {
                if (Double.parseDouble(cutoffFrequency2Field.getText()) <= Double.parseDouble(cutoffFrequency1Field.getText())) {
                    return false;
                }
            }
            if (!testSignalFreqField.getText().isEmpty() && Double.parseDouble(testSignalFreqField.getText()) <= 0) {
                return false;
            }
            if (!testSignalAmpField.getText().isEmpty() && Double.parseDouble(testSignalAmpField.getText()) <= 0) {
                return false;
            }
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    private double parseDouble(String text, double defaultValue) {
        try {
            return text.isEmpty() ? defaultValue : Double.parseDouble(text);
        } catch (NumberFormatException e) {
            return defaultValue;
        }
    }

    private int parseInt(String text, int defaultValue) {
        try {
            return text.isEmpty() ? defaultValue : Integer.parseInt(text);
        } catch (NumberFormatException e) {
            return defaultValue;
        }
    }

    private void showAlert(String message) {
        Alert alert = new Alert(Alert.AlertType.ERROR);
        alert.setTitle("Input Error");
        alert.setHeaderText(null);
        alert.setContentText(message);
        alert.showAndWait();
        LOGGER.warning("Alert shown: " + message);
    }

    private void saveConfig() {
        try {
            FilterConfig config = new FilterConfig(
                    filterType.getValue(), designType.getValue(), orderField.getText(),
                    autoOrderCheck.isSelected(), passbandRippleField.getText(),
                    stopbandAttenuationField.getText(), cutoffFrequency1Field.getText(),
                    cutoffFrequency2Field.getText(), passbandEdgeField.getText(),
                    stopbandEdgeField.getText()
            );
            FileChooser fileChooser = new FileChooser();
            fileChooser.setTitle("Save Filter Configuration");
            fileChooser.getExtensionFilters().add(new FileChooser.ExtensionFilter("JSON Files", "*.json"));
            File file = fileChooser.showSaveDialog(saveButton.getScene().getWindow());
            if (file != null) {
                FilterConfig.saveConfig(config, file);
                showAlert("Configuration saved successfully!");
                statusLabel.setText("Configuration saved");
                LOGGER.info("Configuration saved to: " + file.getAbsolutePath());
            }
        } catch (Exception e) {
            showAlert("Error saving configuration: " + e.getMessage());
            statusLabel.setText("Error saving configuration");
            LOGGER.severe("Error saving configuration: " + e.getMessage());
        }
    }

    private void toggleSweep() {
        if (isSweeping.get()) {
            sweepTimeline.stop();
            isSweeping.set(false);
            startSimButton.setText("Start Sweep");
            statusLabel.setText("Sweep stopped");
            LOGGER.info("Frequency sweep stopped");
        } else {
            isSweeping.set(true);
            startSimButton.setText("Stop Sweep");
            sweepTimeline = new Timeline(new KeyFrame(Duration.millis(100), e -> {
                double value = cutoffFrequency1Slider.getValue() + 0.1;
                if (value > cutoffFrequency1Slider.getMax()) value = cutoffFrequency1Slider.getMin();
                cutoffFrequency1Slider.setValue(value);
                if (cutoffFrequency2Slider.isVisible()) {
                    cutoffFrequency2Slider.setValue(value * 2);
                }
                LOGGER.fine("Sweep updated cutoff to: " + value);
            }));
            sweepTimeline.setCycleCount(Timeline.INDEFINITE);
            sweepTimeline.play();
            statusLabel.setText("Sweeping cutoff frequency...");
            LOGGER.info("Frequency sweep started");
        }
    }

    private void resetFields() {
        filterType.setValue("Type I");
        designType.setValue("Low-pass");
        designMethod.setValue("ZPK");
        autoOrderCheck.setSelected(false);
        orderField.setText("4");
        orderSlider.setValue(4);
        passbandRippleField.setText("1.0");
        passbandRippleSlider.setValue(1.0);
        stopbandAttenuationField.setText("20.0");
        stopbandAttenuationSlider.setValue(20.0);
        cutoffFrequency1Field.setText("1.0");
        cutoffFrequency1Slider.setValue(1.0);
        cutoffFrequency2Field.setText("2.0");
        cutoffFrequency2Slider.setValue(2.0);
        passbandEdgeField.setText("1.0");
        passbandEdgeSlider.setValue(1.0);
        stopbandEdgeField.setText("2.0");
        stopbandEdgeSlider.setValue(2.0);
        testSignalType.setValue("Sine");
        testSignalFreqField.setText("1.0");
        testSignalFreqSlider.setValue(1.0);
        testSignalAmpField.setText("1.0");
        testSignalAmpSlider.setValue(1.0);
        triggerUpdate();
        statusLabel.setText("Fields reset");
        stabilityLabel.setText("Unknown");
        LOGGER.info("Fields reset to default values");
    }

    private void applyTestSignal() {
        if (currentFilter == null) return;
        try {
            double freq = parseDouble(testSignalFreqField.getText(), 1.0);
            double amp = parseDouble(testSignalAmpField.getText(), 1.0);
            String signalType = testSignalType.getValue();
            double duration = 1.0;
            double dt = 0.001;
            int points = (int) (duration / dt) + 1;
            double[] time = generateTimeVector(duration, dt);
            double[] inputSignal = generateTestSignal(signalType, freq, amp, time);
            double[] outputSignal = currentFilter.applyFilter(inputSignal, dt);

            charts.updateSignalChart(time, inputSignal, outputSignal, true);
            statusLabel.setText("Signal applied");
            LOGGER.info("Test signal applied: " + signalType + ", freq=" + freq + ", amp=" + amp);
        } catch (Exception e) {
            showAlert("Error applying test signal: " + e.getMessage());
            statusLabel.setText("Error applying signal");
            LOGGER.severe("Error applying test signal: " + e.getMessage());
        }
    }

    private double[] generateTestSignal(String type, double freq, double amp, double[] time) {
        double[] signal = new double[time.length];
        for (int i = 0; i < time.length; i++) {
            double t = time[i];
            switch (type) {
                case "Sine":
                    signal[i] = amp * Math.sin(2 * Math.PI * freq * t);
                    break;
                case "Square":
                    signal[i] = amp * Math.signum(Math.sin(2 * Math.PI * freq * t));
                    break;
                case "Triangle":
                    signal[i] = amp * (2 * Math.abs(2 * (freq * t - Math.floor(freq * t + 0.5))) - 1);
                    break;
                case "Sawtooth":
                    signal[i] = amp * (2 * (freq * t - Math.floor(freq * t)) - 1);
                    break;
                default:
                    signal[i] = 0;
            }
        }
        return signal;
    }

    private double[] generateTimeVector(double duration, double dt) {
        int points = (int) (duration / dt) + 1;
        double[] time = new double[points];
        for (int i = 0; i < points; i++) {
            time[i] = i * dt;
        }
        return time;
    }
}