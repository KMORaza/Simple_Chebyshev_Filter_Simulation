package simulation.software.codebase;

import com.google.gson.Gson;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class FilterConfig {
    private String filterType;
    private String designType;
    private String order;
    private boolean autoOrder;
    private String passbandRipple;
    private String stopbandAttenuation;
    private String cutoffFrequency1;
    private String cutoffFrequency2;
    private String passbandEdge;
    private String stopbandEdge;

    public FilterConfig(String filterType, String designType, String order, boolean autoOrder,
                        String passbandRipple, String stopbandAttenuation, String cutoffFrequency1,
                        String cutoffFrequency2, String passbandEdge, String stopbandEdge) {
        this.filterType = filterType;
        this.designType = designType;
        this.order = order;
        this.autoOrder = autoOrder;
        this.passbandRipple = passbandRipple;
        this.stopbandAttenuation = stopbandAttenuation;
        this.cutoffFrequency1 = cutoffFrequency1;
        this.cutoffFrequency2 = cutoffFrequency2;
        this.passbandEdge = passbandEdge;
        this.stopbandEdge = stopbandEdge;
    }

    public static void saveConfig(FilterConfig config, File file) throws IOException {
        Gson gson = new Gson();
        try (FileWriter writer = new FileWriter(file)) {
            gson.toJson(config, writer);
        }
    }
}