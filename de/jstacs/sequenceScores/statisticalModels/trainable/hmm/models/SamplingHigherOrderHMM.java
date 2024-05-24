/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models;

import java.io.IOException;
import java.util.Arrays;

import de.jstacs.NotTrainedException;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sampling.BurnInTest;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.SamplingState;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.SimpleSamplingState;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.SamplingEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.SamplingHMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.SamplingTransition;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.TransitionWithSufficientStatistic;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements.TransitionElement;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;

/**
 * 
 * @author Michael Scharfe, Jens Keilwagen
 */
public class SamplingHigherOrderHMM extends HigherOrderHMM {
	
    /**
     * This variable holds the BurnInTest used for training the model
     */
    protected BurnInTest burnInTest;
    /**
     * This boolean indicates if the parameters for the model were sampled
     */
    protected boolean hasSampled;
    private int numberOfStarts;
    private IntList path;
	
    /**
     * This is the main constructor.
     *
     * @param trainingParameterSet the {@link de.jstacs.parameters.ParameterSet} that determines the training algorithm and contains the necessary {@link de.jstacs.parameters.Parameter}s
     * @param name the names of the states
     * @param emissionIdx the indices of the emissions that should be used for each state
     * @param forward a boolean array that indicates whether the symbol on the forward or the reverse complementary strand should be used
     * @param emission the emissions
     * @param te the {@link TransitionElement}s building a transition
     *
     * @throws Exception if
     * 	<ul>
     *  <li>some component could not be cloned</li>
     *  <li>some the length of <code>name, emissionIdx,</code> or <code>forward</code> is not equal to the number of states</li>
     *  <li>not all emissions use the same {@link de.jstacs.data.AlphabetContainer}</li>
     *  <li>the states can not be handled by the transition
     *  </ul>
     */
    public SamplingHigherOrderHMM( SamplingHMMTrainingParameterSet trainingParameterSet, String[] name, int[] emissionIdx, boolean[] forward, SamplingEmission[] emission, TransitionElement... te ) throws Exception {
            super( trainingParameterSet, name, emissionIdx, forward, emission, te );
            hasSampled = false;
            path = new IntList();
    }
	
    /**
     * The standard constructor for the interface {@link de.jstacs.Storable}.
     * Constructs an {@link SamplingHigherOrderHMM} out of an XML representation.
     *
     * @param xml
     *            the XML representation as {@link StringBuffer}
     *
     * @throws NonParsableException
     *             if the {@link SamplingHigherOrderHMM} could not be reconstructed out of
     *             the {@link StringBuffer} <code>xml</code>
     */
    public SamplingHigherOrderHMM( StringBuffer xml ) throws NonParsableException {
            super( xml );
    }
	
    @Override
    public SamplingHigherOrderHMM clone() throws CloneNotSupportedException {
            SamplingHigherOrderHMM clone = (SamplingHigherOrderHMM) super.clone();
            clone.path = path.clone();
            if( burnInTest != null ) {
                    clone.burnInTest = burnInTest.clone();
            }
            return clone;
    }
	
    @Override
    protected void createStates() {
    	this.states = new SimpleSamplingState[emissionIdx.length];
        for( int i = 0; i < emissionIdx.length; i++ ) {
                this.states[i] = new SimpleSamplingState( (SamplingEmission)emission[emissionIdx[i]], name[i], forward[i] );
        }
    }

    /**
     * This method can be used to accept the current parameters (and save them into a file)
     *
     * @throws IOException if the parameters could not be written
     */
    protected void acceptParameters() throws IOException {
    ((SamplingTransition) transition).acceptParameters();
    for( int e = 0; e < emission.length; e++ ) {
            ((SamplingEmission)emission[e]).acceptParameters();
    }
    }
	
    /**
     * This method draws all parameters for the current statistics
     *
     * @throws Exception if the parameters could not be drawn
     */
    protected void drawFromStatistics() throws Exception {
            ((SamplingTransition) transition).drawParametersFromStatistic();
            for( int e = 0; e < emission.length; e++ ) {
                    ((SamplingEmission)emission[e]).drawParametersFromStatistic();
            }
    }

    /**
     * This method implements a sampling step in the sampling procedure
     *
     * @param startPos the start position in the sequence
     * @param endPos the end position in the sequence
     * @param weight the weight for the sequence
     * @param seq the sequence
     * @return the score for the sampi
     * @throws Exception if the sampling step did not succeed
     */
    protected double gibbsSampling( int startPos, int endPos, double weight, Sequence seq ) throws Exception {
            samplePath( path, startPos, endPos, seq );
            addToStatistics( startPos, weight, seq, path );
            return bwdMatrix[0][0];
    }

    private void addToStatistics( int startPos, double weight, Sequence seq, IntList p ) throws Exception{
		int l = 0, layer = 0, state;
		container[1] = 0;

		while( l < p.length() ) {
			state = p.get( l );
			
			int childIdx = transition.getChildIdx( layer, container[1], state );
			if( childIdx < 0 ) {
				throw new IllegalArgumentException( "Impossible path" );
			}
			((SamplingState)states[state]).addToStatistic( startPos, startPos, weight, seq ); //emission
			((SamplingTransition)transition).addToStatistic( layer, container[1], childIdx, weight, seq, startPos ); //transition
			transition.fillTransitionInformation( layer, container[1], childIdx, container );
			if( container[2] == 1 ) {
				startPos++;
				layer++;
			}
			l++;
		}
    }

    private double getLogGammaScoreForCurrentStatistics() {
        double score = ((TransitionWithSufficientStatistic)transition).getLogGammaScoreFromStatistic();
        for(int state = 0; state < states.length; state++) {
            score += ((SamplingState)states[state]).getLogGammaScoreForCurrentStatistic();
        }
        return score;
    }
    
    /**
     * This method implements the next step(s) in the sampling procedure
     * 
     * @param sampling the index of the sampling
     * @param steps the number of sampling that should be executed
     * @param append  whether to append the sampled parameters to an existing file
     *            or to overwrite the file
     * @param data the data used for sampling
     * @param weights the weight for each sequence
     * @throws Exception if something wents wrong
     */
    protected void gibbsSamplingStep( int sampling, int steps, boolean append, DataSet data, double[] weights ) throws Exception {
        double score = 0, weight = 1;
        int N = data.getNumberOfElements();
        Sequence seq;

        //preparation
        burnInTest.setCurrentSamplingIndex(sampling);
        ((SamplingTransition)transition).extendSampling( sampling, append );
        for( int e = 0; e < emission.length; e++ ) {
                ((SamplingEmission)emission[e]).extendSampling( sampling, append );
        }
        sostream.writeln( sampling + " ----------------------------------------" );

        //sampling
        for( int s = 0; s < steps; s++ ) {
                //fill statistics
                resetStatistics();
                score = getLogPriorTerm();
                for( int n = 0; n < N; n++ ) {
                        seq = data.getElementAt( n );
                        if( weights != null ) {
                                weight = weights[n];
                        }
                        score += gibbsSampling( 0, seq.getLength()-1, weight, seq );
                };
                sostream.writeln( s + "\t" + score );
                burnInTest.setValue( score );

                //draw new parameters
                getNewParameters();

                //save the parameters
                acceptParameters();
        }
    }
	
    /**
     * This method set all parameters for the next sampling step
     *
     * @throws Exception if something went wrong
     */
    protected void getNewParameters() throws Exception {
		drawFromStatistics();
    }
	
    @Override
    public void train( DataSet data, double[] weights ) throws Exception {
        int	numberOfBurnInSteps,
                numberOfSteps = ((SamplingHMMTrainingParameterSet)trainingParameter).getNumberOfStepsInStationaryPhase(),
                steps = ((SamplingHMMTrainingParameterSet)trainingParameter).getNumberOfStepsPerIteration(),
                samplingCounter = 0;
        boolean append = false;

        initTraining(data, weights);

        //BURN-IN
        sostream.writeln("GIBBS-SAMPLING - Burn-In ==============================" );
        do {
                for(int start = 0; start < numberOfStarts; start++ ) {
                        if( samplingCounter == 0 ) {
                                initializeRandomly();
                        }

                        gibbsSamplingStep( start, steps, append, data, weights );
                }
                append = true;

                samplingCounter += steps;
                numberOfBurnInSteps = burnInTest.getLengthOfBurnIn();

        } while((samplingCounter-numberOfBurnInSteps) <= 0);

        //FINAL SAMPLING
        sostream.writeln("GIBBS-SAMPLING - Final Sampling =======================");

        for(int start = 0; start < numberOfStarts; start++ ) {
                gibbsSamplingStep( start, numberOfSteps, append, data, weights );
        }

        ((SamplingTransition)transition).samplingStopped();
        for( int e = 0; e < emission.length; e++ ) {
                ((SamplingEmission)emission[e]).samplingStopped();
        }

        hasSampled = true;
    }
	
    @Override
    public String getInstanceName() {
            return "Sampling HMM(" + transition.getMaximalMarkovOrder() + ")";
    }
	
    @Override
    public boolean isInitialized() {
    	return hasSampled;
    }
	
    @Override
    protected double logProb( int startpos, int endpos, Sequence sequence ) throws Exception {
		double logProb = Double.NEGATIVE_INFINITY;
		boolean furtherParam;
		int numSamples = 0;
		
		if( hasSampled ) {

			for(int start = 0; start < numberOfStarts; start++ ) {
				furtherParam = parseParameterSet( start, burnInTest.getLengthOfBurnIn() );
				while( furtherParam ) {
					logProb = Normalisation.getLogSum( logProb, super.logProb( startpos, endpos, sequence ) );
					furtherParam = parseNextParameterSet();
					numSamples++;
				}
			}
		}
		else {
			throw new NotTrainedException();
		}
		
		return logProb - Math.log(numSamples);
	}
	
    /**
     * This method allows the user to parse the set of parameters with index
     * <code>idx</code> of a certain <code>sampling</code> (from a file). The
     * internal numbering should start with 0.
     *
     * @param sampling
     *            the index of the sampling
     * @param idx
     *            the index of the parameter set
     *
     * @return <code>true</code> if the parameter set could be parsed
     *
     * @throws Exception
     *             if there is a problem with parsing the parameters
     */
    protected boolean parseParameterSet( int sampling, int idx ) throws Exception {
            boolean parsed = ((SamplingTransition) transition).parseParameterSet(sampling, idx);
            for( int e = 0; e < emission.length; e++ ) {
                    parsed &= ((SamplingEmission)emission[e]).parseParameterSet(sampling, idx);
            }
            return parsed;
    }
	
    /**
     * This method parse a parameter set stored in file during sampling
     *
     * @return true if parsing was successful
     * @throws Exception if the parameters could not be parsed
     */
    protected boolean parseNextParameterSet() throws Exception {
            boolean parsed = ((SamplingTransition) transition).parseNextParameterSet();
            for( int e = 0; e < emission.length; e++ ) {
                    parsed &= ((SamplingEmission)emission[e]).parseNextParameterSet();
            }
            return parsed;
    }

    private final static String XML_TAG = "SamplingHigherOrderHMM";
	
    @Override
    protected String getXMLTag() {
            return XML_TAG;
    }
	
    @Override
    protected void appendFurtherInformation( StringBuffer xml ) {
            super.appendFurtherInformation( xml );
            XMLParser.appendObjectWithTags( xml, burnInTest, "burnInTest" );
            XMLParser.appendObjectWithTags( xml, hasSampled, "hasSampled" );
    }

    @Override
    protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
            super.extractFurtherInformation( xml );
            numberOfStarts = trainingParameter.getNumberOfStarts();
            burnInTest = XMLParser.extractObjectForTags( xml, "burnInTest", BurnInTest.class );
            hasSampled = XMLParser.extractObjectForTags( xml, "hasSampled", boolean.class );
            path = new IntList();
    }


    /**
     * This methods initialize the training procedure with the given training data
     *
     * @param data the data set used for training
     * @param weights the weight for each sequence
     * @throws Exception if the transition or emissions could not be initialized
     */
    protected void initTraining( DataSet data, double[] weights ) throws Exception {
    	numberOfStarts = trainingParameter.getNumberOfStarts();
   		//GIBBS-SAMPLING INIT
    	burnInTest = ((SamplingHMMTrainingParameterSet)trainingParameter).getBurnInTest();
		burnInTest.resetAllValues();
		((SamplingTransition)transition).initForSampling(numberOfStarts);
		for( int e = 0; e < emission.length; e++ ) {
		        ((SamplingEmission)emission[e]).initForSampling(numberOfStarts);
		}
		furtherInits(data, weights);
    }

    /**
     * This method allows the implementation of further initializations
     *
     * @param data the current data set
     * @param weights the weight for each sequence
     * @throws Exception if the init steps did not succeed
     */
    protected void furtherInits( DataSet data, double[] weights ) throws Exception{
    }

    @Override
    public double[][] getLogStatePosteriorMatrixFor( int startPos, int endPos, Sequence seq ) throws Exception {
        boolean furtherParam;
        double[][] statePosterior = createMatrixForStatePosterior( startPos, endPos );
        double[][] tmp = createMatrixForStatePosterior( startPos, endPos );
        int numSamples = 0;

        for( int s = 0; s < states.length; s++ ) {
            Arrays.fill(statePosterior[s], Double.NEGATIVE_INFINITY);
        }

        if( hasSampled ) {

        	for( int start = 0; start < numberOfStarts; start++ ) {
        		furtherParam = parseParameterSet( start, burnInTest.getLengthOfBurnIn() );
		
                while( furtherParam ) {   
                    fillLogStatePosteriorMatrix( tmp, startPos, endPos, seq, true );
                    for( int s = 0; s < states.length; s++ ) {
                         for( int l = 0; l < statePosterior[s].length; l++ ) {
                             statePosterior[s][l] = Normalisation.getLogSum( statePosterior[s][l], tmp[s][l] );
                         }
                    }
                    furtherParam = parseNextParameterSet();
                    numSamples++;
                }               
            }
        	
        	double d = Math.log( numSamples );
            for( int s = 0; s < states.length; s++ ) {
                 for( int l = 0; l < statePosterior[s].length; l++ ) {
                      statePosterior[s][l] -= d;
                 }
            }
        }    
        else {
        	throw new NotTrainedException();
        }
        
        return getFinalStatePosterioriMatrix( statePosterior );
    }


    @Override
    public double getLogProbForPath( IntList path, int startPos, Sequence seq ) throws Exception {
        boolean furtherParam;
        int numSamples = 0;
        DoubleList d = new DoubleList( 5000 );
        
        if( hasSampled ) {
        	for( int start = 0; start < numberOfStarts; start++ ) {
        		furtherParam = parseParameterSet( start, burnInTest.getLengthOfBurnIn() );
		
                while( furtherParam ) {   
                    d.add( super.getLogProbForPath( path, startPos, seq ) );
                    furtherParam = parseNextParameterSet();
                    numSamples++;
                }               
            }
        	
        	return Normalisation.getLogSum( d.toArray() ) - Math.log( numSamples );
        }    
        else {
        	throw new NotTrainedException();
        }
    }
    
    @Override
    public Pair<IntList,Double> getViterbiPathFor( int startPos, int endPos, Sequence seq ) throws Exception {
        return getViterbiPath( startPos, endPos, seq, ViterbiComputation.MAX );
    }

    /**
     * This method returns a viterbi path that is the optimum for the choosen ViterbiComputation method
     *
     * @param startPos the start position in the sequence
     * @param endPos the end position in the sequence
     * @param seq the sequence
     * @param compute the ViterbiComputation method
     * @return the pair of path and score
     * @throws Exception if the parameters could not be parsed from file
     */
    public Pair<IntList,Double> getViterbiPath( int startPos, int endPos, Sequence seq, ViterbiComputation compute ) throws Exception {
        IntList bestPath = new IntList();
        double bestScore  = Double.NEGATIVE_INFINITY, score;
        
        if(hasSampled) {

            for(int start = 0; start < numberOfStarts; start++ ) {
            	boolean furtherParam = parseParameterSet( start, burnInTest.getLengthOfBurnIn() );
            	
                while( furtherParam ){

                    //SAMPLING VARIANTEN
                    if( compute == ViterbiComputation.SAMPLING ||
                        compute == ViterbiComputation.SAMPLING_GAMMA || 
                        compute == ViterbiComputation.MAX_AND_SAMPLING || 
                        compute == ViterbiComputation.MAX_AND_SAMPLING_GAMMA) {

                    	resetStatistics();
                    	path.clear();

                        gibbsSampling( startPos, endPos, 1d, seq );

                        score = (compute == ViterbiComputation.SAMPLING || compute == ViterbiComputation.MAX_AND_SAMPLING) 
                        	? super.getLogProbForPath( path, startPos, seq )
                        	: getLogGammaScoreForCurrentStatistics();

                        if( score > bestScore) {
                            bestScore = score;
                            bestPath.clear();

                            for( int i =0; i < path.length(); i++ ) {
                                bestPath.add( path.get(i) );
                            }
                        }
                    }

                    //MAXIMIERUNGS VARIANTEN
                    if( compute == ViterbiComputation.MAX ||
                        compute == ViterbiComputation.MAX_GAMMA || 
                        compute == ViterbiComputation.MAX_AND_SAMPLING || 
                        compute == ViterbiComputation.MAX_AND_SAMPLING_GAMMA) {   
                    
                    	resetStatistics();
                    	path.clear();                    	
                        score = viterbi( path, startPos, endPos, 0, seq, null ); //XXX
                        addToStatistics( startPos, 1d, seq, path );

                        score = (compute == ViterbiComputation.MAX || compute == ViterbiComputation.MAX_AND_SAMPLING) 
                        	? score
                        	: getLogGammaScoreForCurrentStatistics();
                        
                        if( score > bestScore) {
                            bestScore = score;
                            bestPath.clear();
                            for( int i =0; i < path.length(); i++ ) {
                            	bestPath.add( path.get(i) );
                            }
                        }
                    }
                    furtherParam = parseNextParameterSet();
                }
            }
        }
        else {
			throw new NotTrainedException();
		}
        return new Pair<IntList, Double>( bestPath, bestScore );
    }

    /**
     * Emumeration of all possible Viterbi-Path methods
     *
     */
    public static enum ViterbiComputation {

        /**
         * standard viterbi path
         */
        MAX,
        /**
         * standard viterbi with gamma scoring (independent from parameters)
         */
        MAX_GAMMA,
        /**
         *  draw paths and score with standard viterbi score
         */
        SAMPLING,
        /**
         * draw paths and score with gamma score (independent from parameters)
         */
        SAMPLING_GAMMA,
        /**
         * best of standard viterbi path and sampled path
         */
        MAX_AND_SAMPLING,
        /**
         *  best of standard viterbi path and sampled path scored with gamma score (independent from parameters)
         */
        MAX_AND_SAMPLING_GAMMA;
    }

    /**
     * This method calculates the a posteriori probability for the current statistics
     *
     * @return the logarithm of the a posteriori probability for the current statistics
     */
    protected double getLogPosteriorFromStatistic() {

        double logPosterior  = ((SamplingTransition)transition).getLogPosteriorFromStatistic();

        for(int state = 0; state < states.length; state++) {
            logPosterior += ((SimpleSamplingState)states[state]).getLogPosteriorFromStatistic();
        }

        return logPosterior;
    }
}