package simulation.software.codebase;

import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.scene.Scene;
import javafx.stage.Stage;
import java.io.IOException;

public class Main extends Application {
    @Override
    public void start(Stage primaryStage) throws IOException {
        FXMLLoader fxmlLoader = new FXMLLoader(Main.class.getResource("/simulation/software/codebase/FilterView.fxml"));
        Scene scene = new Scene(fxmlLoader.load(), 890, 500);
        scene.getStylesheets().add(getClass().getResource("/simulation/software/codebase/styles.css").toExternalForm());
        primaryStage.setTitle("Chebyshev Filter");
        primaryStage.setScene(scene);
        primaryStage.setResizable(false);
        primaryStage.show();
    }

    public static void main(String[] args) {
        launch(args);
    }
}