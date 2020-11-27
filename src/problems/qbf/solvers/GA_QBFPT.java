package problems.qbf.solvers;

import solutions.Solution;

import java.io.IOException;
import java.util.ArrayList;

public class GA_QBFPT extends GA_QBF {

    // array list containing the prohibited triples
    private ArrayList<ArrayList<Integer>> prohibitedTriples = new ArrayList<>();

    /**
     * Constructor for the GA_QBF class. The QBF objective function is passed as
     * argument for the superclass constructor.
     *
     * @param generations  Maximum number of generations.
     * @param popSize      Size of the population.
     * @param mutationRate The mutation rate.
     * @param filename     Name of the file for which the objective function parameters
     *                     should be read.
     * @throws IOException Necessary for I/O operations.
     */
    public GA_QBFPT(Integer generations, Integer popSize, Double mutationRate, String filename) throws IOException {
        super(generations, popSize, mutationRate, filename);

        calculateProhibitedTriples();
    }


    /**
     * Calculate the prohibited triples, using the given functions l(u), g(u) and h(u)
     */
    public void calculateProhibitedTriples() {
        int n = ObjFunction.getDomainSize();

        for (int u = 1; u <= n; u++) {
            ArrayList<Integer> prohibitedTriple = new ArrayList<>();
            prohibitedTriple.add(u);

            // calculate g(u)
            int g_u = calculate_g_u(u, n);

            // calculate h(u)
            int h_u = calculate_h_u(u, n, g_u);

            // add to array
            prohibitedTriple.add(g_u);
            prohibitedTriple.add(h_u);

            // add new prohibited triple to arraylist with all prohibited triples
            prohibitedTriples.add(prohibitedTriple);
        }
    }

    /**
     * Calculates l(u).
     *
     * @param pi_1 first prime number, varies if we are calculating g(u) or h(u)
     * @param pi_2 second prime number, varies if we are calculating g(u) or h(u)
     * @param u    element to which we are finding the other two prohibited values
     * @param n    total number of elements
     * @return result of function l(u)
     */
    public int calculate_l_u(int pi_1, int pi_2, int u, int n) {
        return (1 + ((pi_1 * (u - 1) + pi_2) % n));
    }

    /**
     * Calculates g(u); notice pi_1 and pi_2 values are given on the assignment.
     *
     * @param u element to which we are finding the other two prohibited values
     * @param n total number of elements
     * @return result of function g(u)
     */
    public int calculate_g_u(int u, int n) {
        int pi_1 = 131;
        int pi_2 = 1031;

        int l_u = calculate_l_u(pi_1, pi_2, u, n);

        if (l_u != u)
            return l_u;

        return (1 + (l_u % n));
    }

    /**
     * Calculates h(u); notice pi_1 and pi_2 values are given on the assignment.
     *
     * @param u   element to which we are finding the other two prohibited values
     * @param n   total number of elements
     * @param g_u result of g(u)
     * @return result of function h(u)
     */
    public int calculate_h_u(int u, int n, int g_u) {
        int pi_1 = 193;
        int pi_2 = 1093;

        int l_u = calculate_l_u(pi_1, pi_2, u, n);

        if (l_u != u && l_u != g_u)
            return l_u;
        else if ((1 + (l_u % n)) != u && (1 + (l_u % n)) != g_u) {
            return (1 + (l_u % n));
        }

        return (1 + ((l_u + 1) % n));
    }

    /**
     * We generate a mask with the chromosome size and fill it randomly with 0's and 1's, where a 1
     * means that the alleles are taken from the first parent, while a 0 means they come from the second
     * for each of the two offsprings.
     *
     * @param parents
     * @return
     */
    public Population uniformCrossover(Population parents) {
        Population offsprings = new Population();

        for (int i = 0; i < popSize; i = i + 2) {

            Chromosome parent1 = parents.get(i);
            Chromosome parent2 = parents.get(i + 1);

            ArrayList<Integer> mask = new ArrayList<>();
            for (int j = 0; j < chromosomeSize; j++) {
                mask.add(rng.nextInt(2));
            }

            Chromosome offspring1 = new Chromosome();
            Chromosome offspring2 = new Chromosome();

            for (int j = 0; j < chromosomeSize; j++) {
                if (mask.get(j) == 1) {
                    offspring1.add(parent1.get(j));
                    offspring2.add(parent2.get(j));
                } else {
                    offspring1.add(parent2.get(j));
                    offspring2.add(parent1.get(j));
                }
            }

            offsprings.add(offspring1);
            offsprings.add(offspring2);
        }

        return offsprings;

    }

    /**
     * Generates a random valid chromosome for the QBFPT problem.
     *
     * @return a valid chromosome.
     */
    @Override
    protected Chromosome generateRandomChromosome() {
        Chromosome chromosome = initializeChromosome();
        for (int i = 0; i < chromosomeSize; i++) {
            if (!checkIfValidLocus(chromosome, i))
                chromosome.set(i, 0);
            else
                chromosome.set(i, rng.nextInt(2));
        }

        return chromosome;
    }

    /**
     * Helper function that initializes a chromosome with 0's.
     *
     * @return chromosome with 0's.
     */
    public Chromosome initializeChromosome() {
        Chromosome c = new Chromosome();
        for (int i = 0; i < chromosomeSize; i++) {
            c.add(0);
        }

        return c;
    }

    /**
     * Given a chromosome that may not be valid, i.e. contains prohibited triples, makes it valid
     * by taking off one of the values from each of the existing prohibited triples on the chromosome.
     *
     * @param c
     * @return valid chromosome.
     */
    public Chromosome makeChromosomeValid(Chromosome c) {
        Chromosome validChromosome = c;
        for (int i = 0; i < prohibitedTriples.size(); i++) {
            int firstProhibited = prohibitedTriples.get(i).get(0) - 1;
            int secondProhibited = prohibitedTriples.get(i).get(1) - 1;
            int thirdProhibited = prohibitedTriples.get(i).get(2) - 1;

            // if chromosome has a prohibited triple then take one of the values off
            if (c.get(firstProhibited) == 1 && c.get(secondProhibited) == 1 && c.get(thirdProhibited) == 1) {
                int index = rng.nextInt(3);
                if (index == 0) c.set(firstProhibited, 0);
                if (index == 1) c.set(secondProhibited, 0);
                if (index == 2) c.set(thirdProhibited, 0);
            }
        }

        return validChromosome;
    }

    /**
     * Given a chromosome and a locus, checks if inserting the locus makes the function valid or invalid.
     *
     * @param c
     * @param locus
     * @return true if the function is valid when we add locus; false otherwise.
     */
    public boolean checkIfValidLocus(Chromosome c, Integer locus) {
        for (int i = 0; i < prohibitedTriples.size(); i++) {
            int firstProhibited;
            int secondProhibited;
            int thirdProhibited;
            if (prohibitedTriples.get(i).contains(locus + 1)) {
                /* since we used elements starting with 1 to calculate the prohibited triples,
                 * we need to subtract 1 from each element, since we are starting from zero */
                firstProhibited = prohibitedTriples.get(i).get(0) - 1;
                secondProhibited = prohibitedTriples.get(i).get(1) - 1;
                thirdProhibited = prohibitedTriples.get(i).get(2) - 1;

                if ((firstProhibited == locus && c.get(secondProhibited) == 1 && c.get(thirdProhibited) == 1)
                        || (secondProhibited == locus && c.get(firstProhibited) == 1 && c.get(thirdProhibited) == 1)
                        || (thirdProhibited == locus && c.get(firstProhibited) == 1 && c.get(secondProhibited) == 1))
                    return false;
            }
        }

        return true;
    }

    /**
     * Only mutates the gene if it keeps the chromosome a valid solution to QBFPT.
     *
     * @param chromosome
     * @param locus
     */
    @Override
    protected void mutateGene(Chromosome chromosome, Integer locus) {
        if (chromosome.get(locus) == 1)
            chromosome.set(locus, 0);
        else if (checkIfValidLocus(chromosome, locus))
            chromosome.set(locus, 1);
    }

    /**
     * Method for solving the "default" problem of the assignment.
     *
     * @return
     */
    public Solution<Integer> solveDefault() {

        /* starts the initial population */
        Population population = initializePopulation();

        bestChromosome = getBestChromosome(population);
        bestSol = decode(bestChromosome);
        System.out.println("(Gen. " + 0 + ") BestSol = " + bestSol);

        /*
         * enters the main loop and repeats until a given number of generations
         */
        for (int g = 1; g <= generations; g++) {
            if (checkLimitTime())
                break;

            Population parents = selectParents(population);

            Population offsprings = crossover(parents, popSize);
            for (Chromosome os : offsprings) {
                makeChromosomeValid(os);
            }

            Population mutants = mutate(offsprings);

            Population newpopulation = selectPopulation(mutants);

            population = newpopulation;

            bestChromosome = getBestChromosome(population);

            if (fitness(bestChromosome) > bestSol.cost) {
                bestSol = decode(bestChromosome);
                if (verbose)
                    System.out.println("(Gen. " + g + ") BestSol = " + bestSol);
            }

        }

        printTime();
        return bestSol;
    }

    public Solution<Integer> solveUniformCrossover() {

        /* starts the initial population */
        Population population = initializePopulation();

        bestChromosome = getBestChromosome(population);
        bestSol = decode(bestChromosome);
        System.out.println("(Gen. " + 0 + ") BestSol = " + bestSol);

        /*
         * enters the main loop and repeats until a given number of generations
         */
        for (int g = 1; g <= generations; g++) {
            if (checkLimitTime())
                break;

            Population parents = selectParents(population);

            Population offsprings = uniformCrossover(parents);
            for (Chromosome os : offsprings) {
                makeChromosomeValid(os);
            }

            Population mutants = mutate(offsprings);

            Population newpopulation = selectPopulation(mutants);

            population = newpopulation;

            bestChromosome = getBestChromosome(population);

            if (fitness(bestChromosome) > bestSol.cost) {
                bestSol = decode(bestChromosome);
                if (verbose)
                    System.out.println("(Gen. " + g + ") BestSol = " + bestSol);
            }

        }

        printTime();
        return bestSol;
    }

    /**
     * Randomly select two individuals of the population to be the parents.
     *
     * @param p
     * @return two parents.
     */
    public Population selectParentsSteadyState(Population p) {
        Population parents = new Population();

        int index1 = rng.nextInt(popSize);
        Chromosome parent1 = p.get(index1);
        int index2 = rng.nextInt(popSize);
        Chromosome parent2 = p.get(index2);
        parents.add(parent1);
        parents.add(parent2);

        return parents;
    }

    public Population selectPopulationSteadyState(Population population, Population offsprings) {
        for (int i = 0; i < offsprings.size(); i++) {
            Chromosome worse = getWorseChromosome(population);
            population.remove(worse);
        }
        population.addAll(offsprings);

        return population;
    }

    public Solution<Integer> solveSteadyState() {

        /* starts the initial population */
        Population population = initializePopulation();

        bestChromosome = getBestChromosome(population);
        bestSol = decode(bestChromosome);
        System.out.println("(Gen. " + 0 + ") BestSol = " + bestSol);

        /*
         * enters the main loop and repeats until a given number of generations
         */
        for (int g = 1; g <= generations; g++) {
            if (checkLimitTime())
                break;

            Population parents = selectParentsSteadyState(population);

            Population offsprings = crossover(parents, 2);
            for (Chromosome os : offsprings) {
                makeChromosomeValid(os);
            }

            Population mutants = mutate(offsprings);

            Population newpopulation = selectPopulationSteadyState(population, mutants);

            population = newpopulation;

            bestChromosome = getBestChromosome(population);

            if (fitness(bestChromosome) > bestSol.cost) {
                bestSol = decode(bestChromosome);
                if (verbose)
                    System.out.println("(Gen. " + g + ") BestSol = " + bestSol);
            }

        }

        printTime();
        return bestSol;
    }

    public void printTime() {
        long endTime   = System.currentTimeMillis();
        long totalTime = endTime - startTime;
        System.out.println("Time = "+(double)totalTime/(double)1000+" seg");
    }

    public static void main(String[] args) throws IOException {
        int generations = 10000;
        int popSize1 = 100;
        int popSize2 = 300;
        Double mutationRate1 = 1.0 / 100.0;
        Double mutationRate2 = 10.0 / 100.0;
        String[] instances = {"qbf020", "qbf040", "qbf060", "qbf080", "qbf100", "qbf200", "qbf400"};

        /* PADRÃO: Algoritmo Genético com tamanho de população P1, taxa de mutação M1 e construção
        aleatória da população. */
        System.out.println("1. PADRÃO: Algoritmo Genético com tamanho de população P1, taxa de mutação M1 e construção\n" +
                "aleatória da população.");
        for (int i = 0; i < instances.length; i++) {
            System.out.println("--- INSTANCE " + instances[i] + " ---");
            System.out.println();

            String filename = "instances/" + instances[i];
            GA_QBFPT ga = new GA_QBFPT(generations, popSize1, mutationRate1, filename);
            Solution<Integer> bestSol = ga.solveDefault();
            System.out.println("maxVal = " + bestSol + "\n");
        }

//        /* PADRÃO+POP: Algoritmo Genético PADRÃO mas com tamanho de população P2. */
//        System.out.println("2. PADRÃO+POP: Algoritmo Genético PADRÃO mas com tamanho de população P2.");
//        for (int i = 0; i < instances.length; i++) {
//            System.out.println("--- INSTANCE " + instances[i] + " ---");
//            System.out.println();
//
//            String filename = "instances/" + instances[i];
//            GA_QBFPT ga = new GA_QBFPT(generations, popSize2, mutationRate1, filename);
//            Solution<Integer> bestSol = ga.solveDefault();
//            System.out.println("maxVal = " + bestSol + "\n");
//        }
//
//        /* PADRÃO+MUT: Algoritmo Genético PADRÃO mas com taxa de mutação M2. */
//        System.out.println("3. PADRÃO+MUT: Algoritmo Genético PADRÃO mas com taxa de mutação M2.");
//        for (int i = 0; i < instances.length; i++) {
//            System.out.println("--- INSTANCE " + instances[i] + " ---");
//            System.out.println();
//
//            String filename = "instances/" + instances[i];
//            GA_QBFPT ga = new GA_QBFPT(generations, popSize1, mutationRate2, filename);
//            Solution<Integer> bestSol = ga.solveDefault();
//            System.out.println("maxVal = " + bestSol + "\n");
//        }
//
//        /* PADRÃO+EVOL1: Algoritmo Genético PADRÃO mas com estratégia evolutiva alternativa 1. */
//        System.out.println("4. PADRÃO+EVOL1: Algoritmo Genético PADRÃO mas com estratégia evolutiva alternativa 1.");
//        for (int i = 0; i < instances.length; i++) {
//            System.out.println("--- INSTANCE " + instances[i] + " ---");
//            System.out.println();
//
//            String filename = "instances/" + instances[i];
//            GA_QBFPT ga = new GA_QBFPT(generations, popSize1, mutationRate1, filename);
//            Solution<Integer> bestSol = ga.solveUniformCrossover();
//            System.out.println("maxVal = " + bestSol + "\n");
//        }
//
//        /* PADRÃO+EVOL2: Algoritmo Genético PADRÃO mas com estratégia evolutiva alternativa 2. */
//        System.out.println("5. PADRÃO+EVOL2: Algoritmo Genético PADRÃO mas com estratégia evolutiva alternativa 2.");
//        for (int i = 0; i < instances.length; i++) {
//            System.out.println("--- INSTANCE " + instances[i] + " ---");
//            System.out.println();
//
//            String filename = "instances/" + instances[i];
//            GA_QBFPT ga = new GA_QBFPT(generations, popSize1, mutationRate1, filename);
//            Solution<Integer> bestSol = ga.solveSteadyState();
//            System.out.println("maxVal = " + bestSol + "\n");
//        }

    }

}
