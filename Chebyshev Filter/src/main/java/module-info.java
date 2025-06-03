module simulation.software.codebase.chebyshevfilter {
    requires javafx.controls;
    requires javafx.fxml;
    requires com.google.gson;
    requires java.logging;


    opens simulation.software.codebase to javafx.fxml;
    exports simulation.software.codebase;
}