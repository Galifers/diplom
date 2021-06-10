import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class Main {

    private static XYSeries initialPeak(XYSeries spectrum, int start, int end) {
        XYSeries peak = new XYSeries("Peak");
        for (int j = start; j < end; j++) {
            peak.add(spectrum.getX(j), spectrum.getY(j));
        }
        return peak;
    }

    private static JFreeChart createChart(XYSeriesCollection dataset) {
        JFreeChart chart = ChartFactory
                .createXYLineChart("y = g_1 (x) + ... + g_n (x)", "x", "y",
                        dataset,
                        PlotOrientation.VERTICAL,
                        false, true, true);

        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);

        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setDefaultShapesVisible(false);
        renderer.setDefaultStroke(new BasicStroke(4.0f));
        renderer.setAutoPopulateSeriesStroke(false);

        plot.setRenderer(renderer);
        return chart;
    }

    public static XYSeries loadCSV(String csvFile) {
        String line;
        String cvsSplitBy = " ";
        XYSeries series = new XYSeries("spectrum");
        try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {
            while ((line = br.readLine()) != null) {
                String[] points = line.split(cvsSplitBy);
                series.add(Double.parseDouble(points[0]), Double.parseDouble(points[1]));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return series;
    }

    public static void main(String[] args) {
        XYSeries spectrum = loadCSV("src/main/resources/spectrum.csv");
        XYSeries locations = loadCSV("src/main/resources/locations.csv");
        XYSeriesCollection dataset = new XYSeriesCollection();

        double[] eps = {3.5, 5.5, 3.5, 4.5, 3, 10.4, 22};

        XYSeries peak;

        for (int i = 0; i < locations.getItemCount(); i++) {
            int start = locations.getX(i).intValue();
            int end = locations.getY(i).intValue();
            peak = initialPeak(spectrum, start, end);
            GeneticAlgorithm geneticAlgorithm = new GeneticAlgorithm(8, 1, peak);
            Population population = geneticAlgorithm.initPopulation(3);
            RunnerGA runnerGA = new RunnerGA(geneticAlgorithm, population);
            dataset.addSeries(runnerGA.run(200000, eps[i]));
            peak.clear();
        }

        dataset.addSeries(spectrum);
        JFrame frame = new JFrame("MinimalStaticChart");
        frame.getContentPane().add(new ChartPanel(createChart(dataset)));
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setSize(720, 480);
        frame.setVisible(true);
    }
}