import org.jfree.data.xy.XYSeries;

public class RunnerGA {

    private final GeneticAlgorithm geneticAlgorithm;
    private Population population;

    public RunnerGA(GeneticAlgorithm geneticAlgorithm, Population population) {
        this.geneticAlgorithm = geneticAlgorithm;
        this.population = population;
    }

    private XYSeries deleteLayeringAndBottom(XYSeries series) {
        double maxY = series.getMaxY();
        int maxYi = 0;
        for (int j = 0; j < series.getItemCount(); j++) {
            if ((double) series.getY(j) == maxY) {
                maxYi = j + 1;
                break;
            }
        }
        double sum2 = Double.MAX_VALUE;
        for (int j = maxYi; j < series.getItemCount() - 5; j++) {
            double o1 = (double) series.getY(j);
            double o2 = (double) series.getY(j + 1);
            double o3 = (double) series.getY(j + 2);
            double o4 = (double) series.getY(j + 3);
            double o5 = (double) series.getY(j + 4);
            double sum1 = o1 + o2 + o3 + o4 + o5;
            if (sum1 > sum2) {
                for (int k = j-10; k < series.getItemCount(); k++) {
                    series.update(series.getX(k), -1.0);
                }
                break;
            } else {
                sum2 = sum1;
            }
        }
        int jj = 0;
        int count = -1;
        int itemCount = series.getItemCount();
        while (count++ <= itemCount - 2) {
            if ((double) series.getY(jj) == -1.0 || (double) series.getY(jj) / maxY * 100 < 6.0) {
                series.delete(jj, jj);
                jj--;
            }
            jj++;
        }
        return series;
    }

    private void initialPopulation(Population population, XYSeries series) {
        double maxXSeries = (series.getMinX() + series.getMaxX()) / 2;
        double maxYSeries = series.getMaxY();
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual individual = new Individual(3);
            for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
                double sign = Math.random() >= 0.5 ? -1 : 1;
                double deviation = 0 + Math.random() * 0.01;
                switch (geneIndex) {
                    case 0:
                        individual.setGene(geneIndex, maxYSeries + (sign * (maxYSeries * deviation)));
                        break;
                    case 1:
                        individual.setGene(geneIndex, maxXSeries + (sign * (maxXSeries * deviation)));
                        break;
                    case 2:
                        individual.setGene(geneIndex, 0 + Math.random() * 30);
                        break;
                }
                population.setIndividual(populationIndex, individual);
            }
        }
    }

    public XYSeries run(int plateau, double eps) {

        XYSeries series = geneticAlgorithm.series;

        int start = (int) series.getMinX();
        int end = (int) series.getMaxX();

        series = deleteLayeringAndBottom(series);

        population.initialPopulation(series);

        //initialPopulation(population, series);
        geneticAlgorithm.evalPopulation(population);
        while (!geneticAlgorithm.isTerminationConditionMet(population, eps, plateau)) {
            population = geneticAlgorithm.uniformСrossover(population);
            population = geneticAlgorithm.mutatePopulation(population);
            geneticAlgorithm.evalPopulation(population);
        }
        double a = population.getFittest(0).getGene(0);
        double b = population.getFittest(0).getGene(1);
        double c = population.getFittest(0).getGene(2);
        XYSeries peak = new XYSeries("g(x) = " + a + "*e^-(((x-" + b + "x)^2)/(2*" + c + "))");
        for (int x = start; x < end; x++) {
            double gaussian = a * Math.pow(Math.E, -(Math.pow((x - b), 2.0) / (2 * Math.pow(c, 2.0))));
            if (gaussian > Math.pow(c, 2.2)) {
                peak.add(x, gaussian);
           }
        }
        System.out.println("g(x) = " + a + "*e^-(((x-" + b + ")^2)/(2*" + c + "^2))");
        return peak;
    }
}