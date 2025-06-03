package simulation.software.codebase;

import javafx.animation.FadeTransition;
import javafx.scene.chart.*;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.util.Duration;
import java.util.List;
import java.util.logging.Logger;

public class FilterCharts {
    private static final Logger LOGGER = Logger.getLogger(FilterCharts.class.getName());

    private final TabPane magDbPane, magLinearPane, phasePane, groupDelayPane, poleZeroPane, impulsePane, stepPane, signalPane;
    private final LineChart<Number, Number> magDbChart, magLinearChart, phaseChart, groupDelayChart, impulseChart, stepChart, signalChart;
    private final ScatterChart<Number, Number> poleZeroChart;

    public FilterCharts() {
        // Magnitude (dB)
        NumberAxis magDbX = new NumberAxis(0, 10, 1);
        NumberAxis magDbY = new NumberAxis(-100, 20, 10);
        magDbX.setLabel("Frequency (rad/s)");
        magDbY.setLabel("Gain (dB)");
        magDbChart = new LineChart<>(magDbX, magDbY);
        magDbChart.setCreateSymbols(false);
        magDbPane = new TabPane(new Tab("Magnitude (dB)", magDbChart));

        // Magnitude (Linear)
        NumberAxis magX = new NumberAxis(0, 10, 1);
        NumberAxis magY = new NumberAxis(0, 2, 0.2);
        magX.setLabel("Frequency (rad/s)");
        magY.setLabel("Gain");
        magLinearChart = new LineChart<>(magX, magY);
        magLinearChart.setCreateSymbols(false);
        magLinearPane = new TabPane(new Tab("Magnitude (Linear)", magLinearChart));

        // Phase
        NumberAxis phaseX = new NumberAxis(0, 10, 1);
        NumberAxis phaseY = new NumberAxis(-180, 180, 45);
        phaseX.setLabel("Frequency (rad/s)");
        phaseY.setLabel("Phase (degrees)");
        phaseChart = new LineChart<>(phaseX, phaseY);
        phaseChart.setCreateSymbols(false);
        phasePane = new TabPane(new Tab("Phase", phaseChart));

        // Group Delay
        NumberAxis delayX = new NumberAxis(0, 10, 1);
        NumberAxis delayY = new NumberAxis(0, 100, 10);
        delayX.setLabel("Frequency (rad/s)");
        delayY.setLabel("Group Delay (s)");
        groupDelayChart = new LineChart<>(delayX, delayY);
        groupDelayChart.setCreateSymbols(false);
        groupDelayPane = new TabPane(new Tab("Group Delay", groupDelayChart));

        // Pole-Zero
        NumberAxis pzX = new NumberAxis(-10, 10, 1);
        NumberAxis pzY = new NumberAxis(-10, 10, 1);
        pzX.setLabel("Real");
        pzY.setLabel("Imaginary");
        pzX.setAutoRanging(true);
        pzY.setAutoRanging(true);
        poleZeroChart = new ScatterChart<>(pzX, pzY);
        poleZeroChart.setStyle("-fx-background-color: #ffffff;");
        poleZeroPane = new TabPane(new Tab("Pole-Zero", poleZeroChart));

        // Impulse
        NumberAxis impulseX = new NumberAxis(0, 5, 0.5);
        NumberAxis impulseY = new NumberAxis(-1, 1, 0.2);
        impulseX.setLabel("Time (s)");
        impulseY.setLabel("Amplitude");
        impulseChart = new LineChart<>(impulseX, impulseY);
        impulseChart.setCreateSymbols(false);
        impulsePane = new TabPane(new Tab("Impulse", impulseChart));

        // Step
        NumberAxis stepX = new NumberAxis(0, 5, 0.5);
        NumberAxis stepY = new NumberAxis(-1, 2, 0.5);
        stepX.setLabel("Time (s)");
        stepY.setLabel("Amplitude");
        stepChart = new LineChart<>(stepX, stepY);
        stepChart.setCreateSymbols(false);
        stepPane = new TabPane(new Tab("Step", stepChart));

        // Signal
        NumberAxis signalX = new NumberAxis(0, 1, 0.1);
        NumberAxis signalY = new NumberAxis(-2, 2, 0.5);
        signalX.setLabel("Time (s)");
        signalY.setLabel("Amplitude");
        signalChart = new LineChart<>(signalX, signalY);
        signalChart.setCreateSymbols(false);
        signalPane = new TabPane(new Tab("Signal I/O", signalChart));

        // Enable interactivity
        for (LineChart chart : new LineChart[]{magDbChart, magLinearChart, phaseChart, groupDelayChart, impulseChart, stepChart, signalChart}) {
            chart.setAnimated(true);
            chart.getXAxis().setAutoRanging(false);
            chart.getYAxis().setAutoRanging(true);
        }
        poleZeroChart.setAnimated(true);

        LOGGER.info("Charts initialized");
    }

    public TabPane getMagDbPane() { return magDbPane; }
    public TabPane getMagLinearPane() { return magLinearPane; }
    public TabPane getPhasePane() { return phasePane; }
    public TabPane getGroupDelayPane() { return groupDelayPane; }
    public TabPane getPoleZeroPane() { return poleZeroPane; }
    public TabPane getImpulsePane() { return impulsePane; }
    public TabPane getStepPane() { return stepPane; }
    public TabPane getSignalPane() { return signalPane; }

    public void updateFilterCharts(double[] frequencies,
                                   double[] magDb,
                                   double[] magLinear,
                                   double[] phase,
                                   double[] groupDelay,
                                   List<ChebyshevFilter.Complex> poles,
                                   List<ChebyshevFilter.Complex> zeros,
                                   double[] time,
                                   double[] impulse,
                                   double[] step,
                                   String seriesName,
                                   boolean animate) {
        LOGGER.info("Updating filter charts for series: " + seriesName);

        // Magnitude (dB)
        XYChart.Series<Number, Number> magDbSeries = new XYChart.Series<>();
        magDbSeries.setName(seriesName);
        for (int i = 0; i < frequencies.length && i < magDb.length; i++) {
            if (Double.isFinite(magDb[i])) {
                magDbSeries.getData().add(new XYChart.Data<>(frequencies[i], magDb[i]));
            }
        }
        magDbChart.getData().clear();
        magDbChart.getData().add(magDbSeries);
        if (animate) {
            FadeTransition fade = new FadeTransition(Duration.millis(500), magDbChart);
            fade.setFromValue(0.5);
            fade.setToValue(1.0);
            fade.play();
        }
        LOGGER.fine("Magnitude (dB) chart updated with " + magDbSeries.getData().size() + " points");

        // Magnitude (Linear)
        XYChart.Series<Number, Number> magLinearSeries = new XYChart.Series<>();
        magLinearSeries.setName(seriesName);
        for (int i = 0; i < frequencies.length && i < magLinear.length; i++) {
            if (Double.isFinite(magLinear[i])) {
                magLinearSeries.getData().add(new XYChart.Data<>(frequencies[i], magLinear[i]));
            }
        }
        magLinearChart.getData().clear();
        magLinearChart.getData().add(magLinearSeries);
        if (animate) {
            FadeTransition fade = new FadeTransition(Duration.millis(500), magLinearChart);
            fade.setFromValue(0.5);
            fade.setToValue(1.0);
            fade.play();
        }
        LOGGER.fine("Magnitude (Linear) chart updated with " + magLinearSeries.getData().size() + " points");

        // Phase
        XYChart.Series<Number, Number> phaseSeries = new XYChart.Series<>();
        phaseSeries.setName(seriesName);
        for (int i = 0; i < frequencies.length && i < phase.length; i++) {
            if (Double.isFinite(phase[i])) {
                phaseSeries.getData().add(new XYChart.Data<>(frequencies[i], phase[i]));
            }
        }
        phaseChart.getData().clear();
        phaseChart.getData().add(phaseSeries);
        if (animate) {
            FadeTransition fade = new FadeTransition(Duration.millis(500), phaseChart);
            fade.setFromValue(0.5);
            fade.setToValue(1.0);
            fade.play();
        }
        LOGGER.fine("Phase chart updated with " + phaseSeries.getData().size() + " points");

        // Group Delay
        XYChart.Series<Number, Number> delaySeries = new XYChart.Series<>();
        delaySeries.setName(seriesName);
        for (int i = 0; i < frequencies.length && i < groupDelay.length; i++) {
            if (Double.isFinite(groupDelay[i])) {
                delaySeries.getData().add(new XYChart.Data<>(frequencies[i], groupDelay[i]));
            }
        }
        groupDelayChart.getData().clear();
        groupDelayChart.getData().add(delaySeries);
        if (animate) {
            FadeTransition fade = new FadeTransition(Duration.millis(500), groupDelayChart);
            fade.setFromValue(0.5);
            fade.setToValue(1.0);
            fade.play();
        }
        LOGGER.fine("Group delay chart updated with " + delaySeries.getData().size() + " points");

        // Pole-Zero
        XYChart.Series<Number, Number> stablePoleSeries = new XYChart.Series<>();
        stablePoleSeries.setName("Stable Poles");
        XYChart.Series<Number, Number> unstablePoleSeries = new XYChart.Series<>();
        unstablePoleSeries.setName("Unstable Poles");
        XYChart.Series<Number, Number> zeroSeries = new XYChart.Series<>();
        zeroSeries.setName("Zeros");

        LOGGER.fine("Processing " + poles.size() + " poles and " + zeros.size() + " zeros");
        for (ChebyshevFilter.Complex pole : poles) {
            if (pole != null && Double.isFinite(pole.getReal()) && Double.isFinite(pole.getImag())) {
                XYChart.Series<Number, Number> targetSeries = pole.getReal() < 0 ? stablePoleSeries : unstablePoleSeries;
                targetSeries.getData().add(new XYChart.Data<>(pole.getReal(), pole.getImag()));
                LOGGER.finer("Added pole at (" + pole.getReal() + ", " + pole.getImag() + ") to " + targetSeries.getName());
            }
        }
        for (ChebyshevFilter.Complex zero : zeros) {
            if (zero != null && Double.isFinite(zero.getReal()) && Double.isFinite(zero.getImag())) {
                zeroSeries.getData().add(new XYChart.Data<>(zero.getReal(), zero.getImag()));
                LOGGER.finer("Added zero at (" + zero.getReal() + ", " + zero.getImag() + ")");
            }
        }

        poleZeroChart.getData().clear();
        poleZeroChart.getData().addAll(stablePoleSeries, unstablePoleSeries, zeroSeries);

        // Apply styles
        applySeriesStyle(stablePoleSeries, "-fx-background-color: blue; -fx-shape: \"M0,0 L5,5 L-5,5 L5,-5 L-5,-5 Z\";");
        applySeriesStyle(unstablePoleSeries, "-fx-background-color: red; -fx-shape: \"M0,0 L5,5 L-5,5 L5,-5 L-5,-5 Z\";");
        applySeriesStyle(zeroSeries, "-fx-background-color: green; -fx-shape: \"M-5,-5 L5,5 M-5,5 L5,-5\";");

        if (animate) {
            FadeTransition fade = new FadeTransition(Duration.millis(500), poleZeroChart);
            fade.setFromValue(0.5);
            fade.setToValue(1.0);
            fade.play();
        }
        LOGGER.fine("Pole-zero chart updated with " + stablePoleSeries.getData().size() + " stable poles, " +
                unstablePoleSeries.getData().size() + " unstable poles, and " + zeroSeries.getData().size() + " zeros");

        // Impulse
        XYChart.Series<Number, Number> impulseSeries = new XYChart.Series<>();
        impulseSeries.setName("Impulse");
        for (int i = 0; i < time.length && i < impulse.length; i++) {
            if (Double.isFinite(impulse[i])) {
                impulseSeries.getData().add(new XYChart.Data<>(time[i], impulse[i]));
            }
        }
        impulseChart.getData().clear();
        impulseChart.getData().add(impulseSeries);
        if (animate) {
            FadeTransition fade = new FadeTransition(Duration.millis(500), impulseChart);
            fade.setFromValue(0.5);
            fade.setToValue(1.0);
            fade.play();
        }
        LOGGER.fine("Impulse chart updated with " + impulseSeries.getData().size() + " points");

        // Step
        XYChart.Series<Number, Number> stepSeries = new XYChart.Series<>();
        stepSeries.setName("Step");
        for (int i = 0; i < time.length && i < step.length; i++) {
            if (Double.isFinite(step[i])) {
                stepSeries.getData().add(new XYChart.Data<>(time[i], step[i]));
            }
        }
        stepChart.getData().clear();
        stepChart.getData().add(stepSeries);
        if (animate) {
            FadeTransition fade = new FadeTransition(Duration.millis(500), stepChart);
            fade.setFromValue(0.5);
            fade.setToValue(1.0);
            fade.play();
        }
        LOGGER.fine("Step chart updated with " + stepSeries.getData().size() + " points");

        LOGGER.info("Filter charts updated");
    }

    private void applySeriesStyle(XYChart.Series<Number, Number> series, String style) {
        if (series.getData().isEmpty()) return;
        if (series.getNode() != null) {
            series.getNode().lookupAll(".chart-symbol").forEach(node -> node.setStyle(style));
        }
    }

    public void updateSignalChart(double[] time, double[] inputSignal, double[] outputSignal, boolean animate) {
        XYChart.Series<Number, Number> inputSeries = new XYChart.Series<>();
        inputSeries.setName("Input");
        XYChart.Series<Number, Number> outputSeries = new XYChart.Series<>();
        outputSeries.setName("Output");
        for (int i = 0; i < time.length && i < inputSignal.length && i < outputSignal.length; i++) {
            inputSeries.getData().add(new XYChart.Data<>(time[i], inputSignal[i]));
            outputSeries.getData().add(new XYChart.Data<>(time[i], outputSignal[i]));
        }
        signalChart.getData().clear();
        signalChart.getData().addAll(inputSeries, outputSeries);
        if (animate) {
            FadeTransition fade = new FadeTransition(Duration.millis(500), signalChart);
            fade.setFromValue(0.5);
            fade.setToValue(1.0);
            fade.play();
        }
        LOGGER.info("Signal chart updated with " + inputSeries.getData().size() + " points");
    }
}