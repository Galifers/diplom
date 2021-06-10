import org.jfree.data.xy.XYSeries;

import java.util.Random;

public class GeneticAlgorithm {
    private final int populationSize;

    private double mutationRate = 0.5;

    private final int elitismCount;

    public XYSeries series;

    private int countPlateau = 0;

    private double bestFitness = Double.MAX_VALUE;

    public GeneticAlgorithm(int populationSize, int elitismCount, XYSeries series) {
        this.populationSize = populationSize;
        this.elitismCount = elitismCount;
        this.series = series;
    }

    public Population initPopulation(int chromosomeLength) {
        return new Population(this.populationSize, chromosomeLength);
    }

    public double calcFitness(Individual individual) {
        double fitness = 0;
        double a = individual.getGene(0);
        double b = individual.getGene(1);
        double c = individual.getGene(2);
        for (int i = 0; i < series.getItemCount(); i++) {
            double gaussian = a * Math.pow(Math.E, -(Math.pow(((double) series.getX(i) - b), 2.0) / (2 * Math.pow(c, 2.0))));
            fitness += Math.abs((double) series.getY(i) - gaussian);
        }
        individual.setFitness(fitness);
        return fitness;
    }

    public void evalPopulation(Population population) {
        double populationFitness = 0;
        for (Individual individual : population.getIndividuals()) {
            populationFitness += calcFitness(individual);
        }
        population.setPopulationFitness(populationFitness);
    }

    public double calculateExperientialEps(Population population) {
        double a = population.getFittest(0).getGene(0);
        double b = population.getFittest(0).getGene(1);
        double c = population.getFittest(0).getGene(2);
        double[] raznica = new double[series.getItemCount()];
        double sum = 0;
        for (int i = 0; i < series.getItemCount(); i++) {
            double x = (double) series.getX(i);
            double gaussian = a * Math.pow(Math.E, -(Math.pow((x - b), 2.0) / (2 * Math.pow(c, 2.0))));
            double y = (double) series.getY(i);
            raznica[i] = Math.abs(y - gaussian);
            raznica[i] = raznica[i] / (double) series.getY(i);
            sum += raznica[i];
        }
        return 1.0 / series.getItemCount() * sum * 100;
    }

        public boolean isTerminationConditionMet(Population population, double eps, int plateau) {
            double experientialEps = calculateExperientialEps(population);
            if (experientialEps <= eps || countPlateau == plateau) {
                System.out.printf("eps = %f%n", experientialEps);
                return true;
            } else {
                if (bestFitness == population.getFittest(0).getFitness()) {
                    countPlateau++;
                    mutationRate = (countPlateau / (double) plateau) + 0.5;
                } else {
                    bestFitness = Math.min(bestFitness, population.getFittest(0).getFitness());
                    countPlateau = 0;
                }
            }
            return false;
        }

    private Individual selectParent(Population population) {
        Individual[] individuals = population.getIndividuals();
        double populationFitness = population.getPopulationFitness();
        double[] probabilityOfSelection = new double[individuals.length];
        for (int individualIndex = 0; individualIndex < population.size(); individualIndex++) {
            probabilityOfSelection[individualIndex] = individuals[individualIndex].getFitness() / populationFitness;
        }
        Random rand = new Random();
        double number = rand.nextDouble();
        for (int i = 0; i < probabilityOfSelection.length - 1; i++) {
            if (number <= probabilityOfSelection[i]) {
                return individuals[i];
            }
            probabilityOfSelection[i + 1] += probabilityOfSelection[i];
        }
        return individuals[probabilityOfSelection.length - 1];
    }

    public Population uniformСrossover(Population population) {
        Population newPopulation = new Population(population.size());
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual parent1 = population.getFittest(populationIndex);
            double crossoverRate = 0.95;
            if (crossoverRate > Math.random() && populationIndex >= this.elitismCount) {
                Individual offspring = new Individual(parent1.getChromosomeLength());
                Individual parent2 = selectParent(population);
                double[] v = new double[parent1.getChromosomeLength()];
                for (int i = 0; i < parent1.getChromosomeLength(); i++) {
                    v[i] = Math.random();
                }
                for (int geneIndex = 0; geneIndex < parent1.getChromosomeLength(); geneIndex++) {
                    if (v[geneIndex] <= 0.5) {
                        offspring.setGene(geneIndex, parent1.getGene(geneIndex));
                    } else {
                        offspring.setGene(geneIndex, parent2.getGene(geneIndex));
                    }
                }
                newPopulation.setIndividual(populationIndex, offspring);
            } else {
                newPopulation.setIndividual(populationIndex, parent1);
            }
        }
        return newPopulation;
    }

    public static double getDelta() {
        int m = 20;
        double delta = 0;
        for (int i = 0; i < m; i++) {
            if (Math.random() <= 0.05) {
                delta += Math.pow(2, -i);
            }
        }
        return delta == 0 ? 1 : delta;
    }

    public Population elitePopulation(Population population) {
        Population newPopulation = new Population(this.populationSize);
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual newIndividual = population.getFittest(0);
            newPopulation.setIndividual(populationIndex, newIndividual);
        }
        return newPopulation;
    }

    public Population mutatePopulation(Population population) {
        int m = 25;
        Population newPopulation = new Population(this.populationSize);
        if (mutationRate > 1) {
            newPopulation = elitePopulation(population);
            m *= mutationRate;
        }
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual individual = population.getFittest(populationIndex);
            if (populationIndex >= this.elitismCount) {
                for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
                    if (mutationRate > Math.random()) {
                        double newGene;
                        double sign = Math.random();
                        if (sign <= 0.5) {
                            newGene = Math.abs(individual.getGene(geneIndex) - 0.5 * (individual.getGene(geneIndex) / m)
                                    * getDelta());
                        } else {
                            newGene = Math.abs(individual.getGene(geneIndex) + 0.5 * (individual.getGene(geneIndex) / m)
                                    * getDelta());
                        }

                        double upperBound = series.getMaxY() + series.getMaxY() * 0.02;
                        double lowerBound = series.getMaxY() - series.getMaxY() * 0.02;

                        if (geneIndex == 0) {
                            if (lowerBound <= newGene && newGene <= upperBound) {
                                individual.setGene(geneIndex, newGene);
                            } else {
                                individual.setGene(geneIndex, series.getMaxY());
                            }
                        } else {
                            individual.setGene(geneIndex, newGene);
                        }

                    }
                }
            }
            newPopulation.setIndividual(populationIndex, individual);
        }
        return newPopulation;
    }
}