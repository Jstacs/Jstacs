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
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions;

import java.io.IOException;
import java.text.NumberFormat;
import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class implements a silent emission which is used to create silent states.
 * 
 * @author Jens Keilwagen
 */
public final class SilentEmission implements DifferentiableEmission, SamplingEmission {
	
	/**
	 * The main constructor.
	 */
	public SilentEmission(){}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link SilentEmission} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 */
	public SilentEmission( StringBuffer xml ){
	}
	
	public SilentEmission clone() throws CloneNotSupportedException {
		return (SilentEmission) super.clone();
	}
	
	public StringBuffer toXML() {
		return null;
	}
	
	public double getLogProbAndPartialDerivationFor( boolean forward, int startPos, int endPos, IntList indices, DoubleList partDer,
			Sequence seq ) throws OperationNotSupportedException {
		return 0;
	}
	
	public double getLogProbFor( boolean forward, int startPos, int endPos, Sequence seq ) throws OperationNotSupportedException {
		return 0;
	}

	public void joinStatistics(Emission... emissions){}
	
	public void addToStatistic( boolean forward, int startPos, int endPos, double weight, Sequence seq ) throws OperationNotSupportedException {}

	public void estimateFromStatistic() {}
	
	public void resetStatistic() {}
	
	public AlphabetContainer getAlphabetContainer() {
		return null;
	}

	public void addGradientOfLogPriorTerm( double[] gradient, int offset ) {}

	public double getLogPriorTerm() {
		return 0;
	}

	public void fillCurrentParameter( double[] params ) {}
	
	public void setParameter( double[] params, int offset ) {}

	public int setParameterOffset( int offset ) {
		return offset;
	}
	
	public void initializeFunctionRandomly() {}

	public void drawParametersFromStatistic() {}

	public void extendSampling( int sampling, boolean append ) throws IOException {}

	public void initForSampling( int starts ) throws IOException {}

	public boolean isInSamplingMode() {
		return true;
	}

	public boolean parseNextParameterSet() {
		return true;
	}

	public boolean parseParameterSet( int sampling, int n ) throws Exception {
		return true;
	}

	public void samplingStopped() throws IOException {}

	public void acceptParameters() throws IOException {}
	
	public double getLogPosteriorFromStatistic() { return 0; }
	
	public String toString( NumberFormat nf ) {
		return this.getClass().getSimpleName() + "\n";
	}

	public double getLogGammaScoreFromStatistic() {
		return 0;
	}

	@Override
	public String getNodeShape(boolean forward) {
		return "\"circle\"";
	}

	@Override
	public String getNodeLabel( double weight, String name, NumberFormat nf ) {
		return "\""+name+"\"";
	}

	@Override
	public void fillSamplingGroups( int parameterOffset, LinkedList<int[]> list ) {
		
	}

	@Override
	public int getNumberOfParameters() {
		return 0;
	}

	@Override
	public int getSizeOfEventSpace() {
		return 0;
	}
	
	@Override
	public void setParameters(Emission t) throws IllegalArgumentException {
		if( !t.getClass().equals( getClass() ) ) {
			throw new IllegalArgumentException( "The emissions are not comparable." );
		}		
	}

	@Override
	public boolean isNormalized() {
		return true;
	}	
}
