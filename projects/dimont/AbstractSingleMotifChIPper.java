/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package projects.dimont;

import java.text.NumberFormat;
import java.util.HashMap;
import java.util.Iterator;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.DataSet;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.NonParsableException;
import de.jstacs.motifDiscovery.Mutable;
import de.jstacs.motifDiscovery.MutableMotifDiscoverer;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.NormalizedDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.VariableLengthDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.AbstractMixtureDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.StrandDiffSM;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.ToolBox;


/**
 * Prototype
 * 
 * implemented
 * <ul>
 * <li>profile</li>
 * <li>thresholded variant</li>
 * <li>{@link MutableMotifDiscoverer}</li>
 * </ul>
 * 
 * not implemented
 * <ul>
 * <li>flanking</li>
 * <li>multiple motif</li>
 * </ul>
 * 
 *   
 * @author Jens Keilwagen
 */
public abstract class AbstractSingleMotifChIPper extends AbstractMixtureDiffSM implements MutableMotifDiscoverer, VariableLengthDiffSM {
	
	protected double logP;
	protected HashMap<Sequence,float[]> positionHash;
	
	public AbstractSingleMotifChIPper( int starts, DifferentiableStatisticalModel motif ) throws CloneNotSupportedException {
		super( 0, starts, 2, true, true, motif );
		init();
	}
	
	public AbstractSingleMotifChIPper( StringBuffer xml ) throws NonParsableException {
		super( xml );
		init();
	}
	
	protected void init() {
		logP = -Math.log( function[0].getAlphabetContainer().getAlphabetLengthAt(0) );
		positionHash = new HashMap<Sequence, float[]>();
	}
	
	public AbstractSingleMotifChIPper clone() throws CloneNotSupportedException {
		AbstractSingleMotifChIPper clone = (AbstractSingleMotifChIPper) super.clone();
		//clone.end = end.clone();
		clone.positionHash = new HashMap<Sequence, float[]>();
		Iterator<Sequence> it = positionHash.keySet().iterator();
		while( it.hasNext() ) {
			Sequence s = it.next();
			clone.positionHash.put( s, positionHash.get(s).clone() );
		}
		return clone;
	}

	protected boolean determineIsNormalized(){
		return false;
	}
	
	protected Sequence getReference( Sequence seq ) {
		SequenceAnnotation seqAn = seq.getSequenceAnnotationByTypeAndIdentifier("reference", "reads");
		return seqAn == null ? null : ((ReferenceSequenceAnnotation) seqAn).getReferenceSequence();
	}
	
	protected float[] getPosition( Sequence seq, boolean add ) {
		float[] res = positionHash.get( seq );
		if( res == null ) {
			Sequence ref = getReference( seq );
			if( ref != null ) {
				res = new float[seq.getLength()-function[0].getLength()+1];
				float sum = 0;
				for( int i = 0; i < res.length; i++ ) {
					res[i] = (float) ref.continuousVal( i );
					sum+= res[i];
				}				
				for( int i = 0; i < res.length; i++ ) {
					res[i] = (float)Math.log( res[i]/sum );
				}				
			}
			if( add ) {
				positionHash.put( seq, res );
			}
		}
		return res;
	}
	
	
	
	
	
	@Override
	public double getLogNormalizationConstant(int length) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getLogPartialNormalizationConstant(int parameterIndex, int length) throws Exception {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setStatisticForHyperparameters(int[] length, double[] weight) throws Exception {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double getLogScoreFor(Sequence seq, int startpos, int endpos){
		return getLogScoreFor( seq.getSubSequence( startpos, endpos-startpos+1 ), 0 );
	}

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int startpos, int endpos, IntList indices,
			DoubleList partialDer) {
		return getLogScoreAndPartialDerivation(seq.getSubSequence( startpos, endpos-startpos+1 ), 0, indices, partialDer);
	}

	@Override
	public double getHyperparameterForHiddenParameter( int index ) {
		//XXX
		return function[0].getESS();
	}
	
	@Override
	protected double getLogNormalizationConstantForComponent(int i) {
		if( i == 0 )
		{
			return function[i].getLogNormalizationConstant();
		}
		else
		{
			return 0;
		}
	}

	public static int draw( DataSet d, double[] weight ) {
		if( weight == null ) {
			return r.nextInt( d.getNumberOfElements() );
		} else {
			return draw( weight );
		}
	}
	
	public static int draw( double[] weight ) {
		double s = r.nextDouble() * ToolBox.sum( weight );
		int i = 0;
		while( weight[i] < s ) {
			s -=weight[i];
			i++;
		}
		//System.out.println( weight[i] + "\t[" + ToolBox.min(weight) + "," + ToolBox.max(weight) + "]" );
		return i;
	}
	
	@Override
	protected void initializeUsingPlugIn(int index, boolean freeParams, DataSet[] data, double[][] weights) throws Exception {
		// consensus 90%
		int a = (int) alphabets.getAlphabetLengthAt(0), s, p, l;
		double d = 0.1/(a-1), h;
		d = (1d-a*d)/(a*d);
		Sequence seq, ref;

		l = function[0].getLength();
		s = draw( data[index], weights == null ? null : weights[index] );
		seq = data[index].getElementAt( s );
		
		ref = getReference(seq);
		if( ref == null ) {
			p = r.nextInt( seq.getLength() - l + 1 );
		} else {
			double[] prof = new double[seq.getLength()-function[0].getLength()];
			for( int i = 0; i < prof.length; i++ ) {
				prof[i] = ref.continuousVal( i );
			}
			p = draw( prof );
		}
		seq = seq.getSubSequence( p, l );
		h = d * function[0].getESS();
		
		this.freeParams = freeParams;
		initializeMotif(index, new DataSet( "", seq ), new double[]{h} );//motifIndex is equal to index. Does not make sense. TODO FIXME
	}
	
	@Override
	public void initializeMotif(int motifIndex, DataSet data, double[] weights) throws Exception {
		int a = data.getElementLength()-function[0].getLength();
		//System.out.println("before: "+function[0].getNumberOfParameters());
		//System.out.println(a);
		if( a != 0 ) {
			((Mutable)function[0]).modify( 0, a );
		}
		//System.out.println("between: "+function[0].getNumberOfParameters());
		function[motifIndex].initializeFunction( 0, freeParams, new DataSet[]{data}, new double[][]{weights} );
		if( a != 0 ) {
			double d = a/2d;
			((Mutable)function[0]).modify( (int)Math.ceil(d), (int)Math.ceil(-d) );
		}
		//System.out.println("after: "+function[0].getNumberOfParameters());
		init( freeParams );
		initializeHiddenUniformly();
		positionHash.clear();
		//System.out.println("finally: "+function[0].getNumberOfParameters());
	}

	@Override
	public double getESS() {
		double ess = 0;
		for( int i = 0; i <= function.length; i++ ) {
			ess += getHyperparameterForHiddenParameter(i);
		}
		return ess;
	}

	@Override
	public double getLogPartialNormalizationConstant(int parameterIndex)
			throws Exception {
		if( Double.isNaN( norm ) )
		{
			precomputeNorm();
		}
		int[] ind = getIndices( parameterIndex );
		double res;
		if( ind[0] == function.length )
		{
			res = partNorm[ind[1]];
		}
		else {
			res = logHiddenPotential[ind[0]] + function[ind[0]].getLogPartialNormalizationConstant( ind[1] );
		}
		
		//System.out.println( parameterIndex + "\t" + Arrays.toString(ind) + "\t" + res );
		
		return res;
	}

	@Override
	public String getInstanceName() {
		return getClass().getSimpleName() + "(" + function[0].getInstanceName() + ")";
	}

	@Override
	public int getGlobalIndexOfMotifInComponent( int component, int motif ) {
		return component;
	}
	
	@Override
	public int getIndexOfMaximalComponentFor(Sequence sequence) throws Exception {
		return getIndexOfMaximalComponentFor( sequence, 0 );
	}

	@Override
	public int getMotifLength(int motif) {
		return function[motif].getLength();
	}

	@Override
	public int getNumberOfMotifs() {
		return function.length;
	}

	@Override
	public int getNumberOfMotifsInComponent(int component) {
		if( component < function.length )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	
	protected abstract int fillMotifComponentScoreOf( double[] array, Sequence sequence, int startpos );

	public abstract double[] getStrandProfileOfScoresFor( Sequence seq, boolean forwardStrand ) throws OperationNotSupportedException;
	
	@Override
	public double[] getProfileOfScoresFor( int component, int motif, Sequence sequence, int startpos, KindOfProfile kind ) throws Exception {
		if( motif == 0 && component < function.length )
		{
			//unnormalized log score
			double d = 0;
			int l = sequence.getLength()- startpos - function[component].getLength();
			double[] res;
			if( l >= 0 ) {
				res = new double[l+1];
				int end = fillMotifComponentScoreOf( res, sequence, startpos );
				switch( kind )
				{
					case UNNORMALIZED_JOINT:
						d = logHiddenPotential[component];
					case UNNORMALIZED_CONDITIONAL:
						//XXX d += l * logP;
						break;				
					case NORMALIZED_CONDITIONAL:
						d = -Normalisation.getLogSum( 0, end, res );
						break;
					
					default:
						throw new IndexOutOfBoundsException();
				}
				for( int i = 0; i < res.length; i++ )
				{
					res[i] += d;
				}
			} else {
				res = new double[0];
			}
			
			return res;
		}
		else
		{
			throw new IndexOutOfBoundsException();
		}
	}

	@Override
	public double[] getStrandProbabilitiesFor(int component, int motif, Sequence sequence, int startpos) throws Exception {
		if( motif > 0 || component > function.length )
		{
			throw new IndexOutOfBoundsException();
		}
		else
		{
			DifferentiableSequenceScore m = function[component];
			while( m instanceof NormalizedDiffSM )
			{
				m = ((NormalizedDiffSM)m).getFunction();
			}
			if( m instanceof StrandDiffSM )
			{
				if(startpos == 0){
					return ((StrandDiffSM)m).getProbsForComponent( sequence );
				}else{
					return ((StrandDiffSM)m).getProbsForComponent( sequence.getSubSequence( startpos ) );
				}
			}
			else
			{
			    return new double[]{1.0,0.0};
			}
		}
	}
	
	
	@Override
	public void adjustHiddenParameters(int index, DataSet[] data, double[][] weights) throws Exception {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void initializeMotifRandomly(int motif) throws Exception {
		function[motif].initializeFunctionRandomly( freeParams );
		init( freeParams );
		
	}
	
	protected void init( boolean freeParams ) {
		super.init( freeParams );
		if( positionHash != null ) {
			positionHash.clear();
		}
	}

	@Override
	public boolean modifyMotif(int motifIndex, int offsetLeft, int offsetRight) throws Exception {
		if( function[motifIndex] instanceof Mutable ) {
			double norm_old = function[motifIndex].getLogNormalizationConstant();
			boolean res = ((Mutable) function[motifIndex]).modify(offsetLeft, offsetRight);
			if( res )
			{
				init( freeParams );
				double norm_new = function[motifIndex].getLogNormalizationConstant();
				hiddenParameter[motifIndex] += ( norm_old - norm_new );
				this.setHiddenParameters( hiddenParameter , 0 );
				norm = Double.NaN;
			}
			
			return res;
		} else {
			return false;
		}
	}
	
	public String toString(NumberFormat nf)
	{
		if( Double.isNaN( norm ) )
		{
			precomputeNorm();
		}
		StringBuffer erg = new StringBuffer( function.length * 1000 );
		erg.append( "\nno motif: " + Math.exp(partNorm[function.length] - norm) + "\texp(" +partNorm[function.length] + " - " + norm + ")\t" + logHiddenPotential[function.length] + "\n" );
		for( int i = 0; i < function.length; i++ ) {
			erg.append( "\nmotif " + i + ": " );
			erg.append( Math.exp(partNorm[i] - norm) +"\texp(" +partNorm[i] + " - " + norm + ")\t" + logHiddenPotential[i] );
			erg.append( "\n" + function[i].toString() + "\n" );
		}
		return erg.toString();
	}
	
	public abstract void reset();
	
	public void resetPositions(){
		positionHash.clear();
	}
}
