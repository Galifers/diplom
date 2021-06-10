public class Individual {

    private double[] chromosome;
    private double fitness = Double.MAX_VALUE;

    public Individual(double[] chromosome) {
        this.chromosome = chromosome;
    }

    public Individual(int chromosomeLength) {
        this.chromosome = new double[chromosomeLength];
        for (int gene = 0; gene < chromosomeLength; gene++) {
            this.setGene(gene, 1 + Math.random() * 10);
        }
    }

    public double[] getChromosome() {
        return this.chromosome;
    }

    public int getChromosomeLength() {
        return this.chromosome.length;
    }

    public void setGene(int offset, double gene) {
        this.chromosome[offset] = gene;
    }

    public double getGene(int offset) {
        return this.chromosome[offset];
    }

    public void setFitness(double fitness) {
        this.fitness = fitness;
    }

    public double getFitness() {
        return this.fitness;
    }

    public String toString() {
        String output = "[";
        for (int gene = 0; gene < this.chromosome.length; gene++) {
            output += this.chromosome[gene] + " ";
        }
        output += "]";
        return output;
    }
}
