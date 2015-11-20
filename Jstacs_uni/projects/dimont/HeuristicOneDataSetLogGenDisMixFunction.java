package projects.dimont;

import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.OneDataSetLogGenDisMixFunction;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior;
import de.jstacs.data.DataSet;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;

public class HeuristicOneDataSetLogGenDisMixFunction extends OneDataSetLogGenDisMixFunction {
	
	public HeuristicOneDataSetLogGenDisMixFunction(int threads, DifferentiableSequenceScore[] score, DataSet data, double[][] weights, LogPrior prior, double[] beta, boolean norm, boolean freeParams) throws IllegalArgumentException {
		super(threads, score, data, weights, prior, beta, norm, freeParams);
	}

	public void resetHeuristics() {
		for( int i = 0; i < score.length; i++ ) {
			for( int j = 0; j < score[i].length; j++ ) {
				if( score[i][j] instanceof AbstractSingleMotifChIPper ) {
					//System.out.println("reset " + i + " " + j );
					((AbstractSingleMotifChIPper) score[i][j]).reset();
				}
			}
		}
	}	
}
