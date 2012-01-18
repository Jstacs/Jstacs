package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models;

import java.util.Random;

import de.jstacs.data.WrongAlphabetException;
import de.jstacs.io.NonParsableException;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.discrete.PhyloDiscreteEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.SamplingHMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.SamplingTransition;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements.TransitionElement;

/**
 * This class implements an (higher order) HMM that contains multi-dimensional emissions described
 * by a phylogenetic tree. The model is trained by a Metropolis-Hastings algorithm.
 *
 * @author Michael Scharfe
 */
public class SamplingPhyloHMM extends SamplingHigherOrderHMM {

    static Random r = new Random();
    
    /**
     * This is the main constructor for a hidden markov model with phylogenetic emission(s)
     * This model can be trained using a metropolis hastings algorithm
     *
     * @param trainingParameterSet the {@link de.jstacs.parameters.ParameterSet} that determines the training algorithm and contains the necessary {@link de.jstacs.parameters.Parameter}s
     * @param name the names of the states
     * @param emissionIdx the indices of the emissions that should be used for each state
     * @param forward a boolean array that indicates whether the symbol on the forward or the reverse complementary strand should be used
     * @param emission the emissions
     * @param te the {@link TransitionElement}s building a transition
     * 
     * @throws CloneNotSupportedException if the parameters or the emissions or transition elements could not be cloned
     * @throws IllegalArgumentException if one of the parameters is not allowed
     * @throws WrongAlphabetException if the alphabet does not fit the model
     * @throws Exception if something else went wrong
     */
    public SamplingPhyloHMM( SamplingHMMTrainingParameterSet trainingParameterSet, String[] name, int[] emissionIdx, boolean[] forward, PhyloDiscreteEmission[] emission, TransitionElement... te  ) throws CloneNotSupportedException, IllegalArgumentException, WrongAlphabetException, Exception {

        super(trainingParameterSet, name, emissionIdx, forward, emission, te);
    }

    /**
     * The standard constructor for the interface {@link de.jstacs.Storable}.
     * Constructs an {@link SamplingPhyloHMM} out of an XML representation.
     *
     * @param xml
     *            the XML representation as {@link StringBuffer}
     *
     * @throws NonParsableException
     *             if the {@link SamplingPhyloHMM} could not be reconstructed out of
     *             the {@link StringBuffer} <code>xml</code>
     */
    public SamplingPhyloHMM( StringBuffer xml ) throws NonParsableException {
            super( xml );
    }

    @Override
    public String getInstanceName() {
        return "PhyloHMM(" + transition.getMaximalMarkovOrder() + ")";
    }

    @Override
    protected void getNewParameters() throws Exception {

        boolean accepted = false;   
        double acceptProb;

        while(!accepted) {

            acceptProb = 0d;
            //current parameters
            acceptProb =  getLogProposalPosteriorFromStatistic() - getLogPosteriorFromStatistic();

            //draw proposal from non-phylo distributions
            drawFromStatistics();

            //proposed parameters
            acceptProb += getLogPosteriorFromStatistic() - getLogProposalPosteriorFromStatistic();
            acceptProb = Math.exp(acceptProb);

            if(acceptProb >= 1d || drawIndexFrom(new double[]{1d-acceptProb, acceptProb}) > 0) {
                acceptParameters();
                accepted = true;
            }
        }
    }

    private double getLogProposalPosteriorFromStatistic() {

        double logPosterior  = ((SamplingTransition)transition).getLogPosteriorFromStatistic();

        for( int e = 0; e < emission.length; e++ ) {
            logPosterior += ((PhyloDiscreteEmission)emission[e]).getLogProposalPosteriorFromStatistic();
	}
        return logPosterior;
    }



    private static int drawIndexFrom(double[] distribution) {

        int index = 0;

        double p = r.nextDouble();

        while(index < distribution.length && p > distribution[index]) {
                p -= distribution[index++];
        }
        if(index == distribution.length ) {
                index--;
        }

        return index;
    }
}
