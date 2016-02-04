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
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.motifDiscovery.MotifDiscoverer.KindOfProfile;
import de.jstacs.motifDiscovery.MutableMotifDiscoverer;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.NormalizedDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSM;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;


/**
 * Prototype
 * 
 * implemented
 * <ul>
 * <li>profile</li>
 * <li>thresholded variant</li>
 * <li>{@link MutableMotifDiscoverer}</li>
 * <li>multiple motif</li>
 * </ul>
 * 
 * not implemented
 * <ul>
 * <li>flanking</li>
 * </ul>
 * 
 *   
 * @author Jens Keilwagen
 */
public class ThresholdedStrandChIPper extends AbstractSingleMotifChIPper {
	
	private static final double LOG2 = Math.log( 2 );
	
	private double[] scoreProfile;
	private double[] normedProfile;
	private double[] temp;
	private IntList tempPos;
	private IntList end;
	private double threshold;
	private HashMap<Sequence,IntList> toBeUsedHash;
	
	public ThresholdedStrandChIPper( int starts, double t, DifferentiableStatisticalModel motif ) throws CloneNotSupportedException {
		super( starts, motif );
		if( t <= 0 || t > 1 ) {
			throw new IllegalArgumentException();
		}
		threshold = t;
		init();
	}
	
	public ThresholdedStrandChIPper( StringBuffer xml ) throws NonParsableException {
		super( xml );
		init();
	}
	
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
		threshold = (Double) XMLParser.extractObjectForTags( xml, "threshold" );
	}
	
	protected StringBuffer getFurtherInformation() {
		StringBuffer extra = new StringBuffer(100);
		XMLParser.appendObjectWithTags( extra, threshold, "threshold" );
		return extra;
	}
	
	protected void init() {
		super.init();
		end = new IntList();
		tempPos = new IntList();
		toBeUsedHash = new HashMap<Sequence, IntList>();
	}
	
	public ThresholdedStrandChIPper clone() throws CloneNotSupportedException {
		ThresholdedStrandChIPper clone = (ThresholdedStrandChIPper) super.clone();
		if(scoreProfile != null){
			clone.scoreProfile = scoreProfile.clone();
		}
		if(normedProfile != null){
			clone.normedProfile = normedProfile.clone();
		}
		if(temp != null){
			clone.temp = temp.clone();
		}
		clone.end = end.clone();
		clone.tempPos = tempPos.clone();
		clone.toBeUsedHash = new HashMap<Sequence, IntList>();
		Iterator<Sequence> it = toBeUsedHash.keySet().iterator();
		while( it.hasNext() ) {
			Sequence s = it.next();
			clone.toBeUsedHash.put( s, toBeUsedHash.get( s ).clone() );
		}
		return clone;
	}
	
	
	protected int fillMotifComponentScoreOf( double[] array, Sequence sequence, int startpos ) {
		int end = fillComponentScoreOf( array, sequence, startpos, null, false );
		return end;
	}
	
	protected int fillComponentScoreOf( Sequence seq, int start, IntList b, boolean usePosition )
	{
		int m = function[0].getLength();
		int end = seq.getLength()-start-m+1;
		if( scoreProfile == null || scoreProfile.length < 2*end ) {
			scoreProfile = new double[seq.getLength()*2];
		}
		return fillComponentScoreOf( scoreProfile, seq, start, b, usePosition );
	}
	
	protected int fillComponentScoreOf( double[] array, Sequence seq, int start, IntList b, boolean usePosition ) {
		//System.out.println(seq);
		int m = function[0].getLength();
		int end = seq.getLength()-start-m+1;
		double pos = 0;
		float[] position = null;
		if( usePosition ) {
			position = getPosition( seq, false );
			if( position == null ) {
				pos = - Math.log(end);
			}
		}
		int j = 0, curr = -1, n = b == null ? end : b.length();
		for(int l=0;l<n;l++){
			if(b != null){
				curr = b.get( l );
				boolean rc = false;
				if(curr < 0){
					rc = true;
					curr = -curr-1;
				}
				if( curr >= start && curr <= seq.getLength()-m ){
					if( position != null ) {
						pos = position[curr];//XXX
					}
					//System.out.print(pos+" ");
					if(rc){
						try {
							array[j++] = pos -LOG2 + function[0].getLogScoreFor( seq.reverseComplement(), seq.getLength()-curr-m ) - m*logP;
						} catch ( OperationNotSupportedException e ) {
							throw new RuntimeException( e );
						}
					}else{
						array[j++] = pos -LOG2 + function[0].getLogScoreFor( seq, curr ) - m*logP;
					}
				}
				
			}else{
				curr = start+l;
				if( curr >= start && curr <= seq.getLength()-m ){
					if( position != null ) {
						pos = position[curr];//XXX
					}
					//System.out.print(pos+" ");
					try {
						array[j++] = Normalisation.getLogSum( pos -LOG2 + function[0].getLogScoreFor( seq, curr ) - m*logP,
								pos -LOG2 + function[0].getLogScoreFor( seq.reverseComplement(), seq.getLength()-curr-m ) - m*logP);
					} catch ( OperationNotSupportedException e ) {
						throw new RuntimeException( e );
					}
					
				}
			}
			
		}
		//System.out.println();
		return j;
	}
	
	public double[] getStrandProfileOfScoresFor( Sequence seq, boolean forwardStrand ) throws OperationNotSupportedException {
		int m = function[0].getLength();
		int end = seq.getLength()-m+1;
		double[] res = new double[end];
		double offset = logHiddenPotential[0] - Math.log(end) + end*logP;
		if( forwardStrand ) {
			for ( int i = 0; i < end; i++ ) {
				res[i] = offset + function[0].getLogScoreFor( seq, i );
			}
		} else {
			for ( int i = 0; i < end; i++ ) {
				res[i] = offset + function[0].getLogScoreFor( seq.reverseComplement(), end-i-1 );
			}
		}
		
		return res;
	}

	protected void fillComponentScores( Sequence seq, int start )
	{
		double h = (seq.getLength() - start)*logP;
		IntList b = toBeUsedHash.get(seq);
		int j = fillComponentScoreOf( seq, start, b==null ? null : b, true );
		componentScore[0] = h + logHiddenPotential[0] + Normalisation.getLogSum( 0, j, scoreProfile );
		componentScore[function.length] = h + logHiddenPotential[function.length];
	}

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList partialDer) {
		int l, m, n, counter, stop;
		int startIndex = partialDer.length();
		
		double pos = Double.NaN;
		
		IntList b = toBeUsedHash.get(seq);
		boolean add = b == null;
		if( add ) {
			b = new IntList();
		} else {
			//System.out.println("size "+b.length());
			b.clear();
		}
		
		// for each motif (function)
		m = function[0].getLength();
		stop = seq.getLength() - start - m + 1;
		if( stop > 0 ) {
			float[] position = getPosition( seq, true );
			if( position == null ) {
				pos = - Math.log(stop);
			}

			if( scoreProfile == null || scoreProfile.length < 2*stop ) {
				scoreProfile = new double[seq.getLength()*2];
				
			}
			if( normedProfile == null || normedProfile.length < 2*stop ) {
				normedProfile = new double[seq.getLength()*2];
			}
			if( temp == null || temp.length < 2*stop ) {
				temp = new double[seq.getLength()*2];
			}
			
			
			// for each start position
			for( n = 0, l = start; n < stop; l++, n++ ) {
				if( position != null ) {
					pos = position[l];//XXX
				}
				
				// compute the score
				scoreProfile[n] =
					pos // position
					- LOG2 + function[0].getLogScoreFor( seq, l ) // motif
					- m*logP
					;
				try {
					scoreProfile[n + stop] = 
						pos // position
						- LOG2 + function[0].getLogScoreFor( seq.reverseComplement(), seq.getLength() - l - m ) // motif
						- m*logP
						;
				} catch ( OperationNotSupportedException e1 ) {
					throw new RuntimeException( e1 );
				}
			}
			
			Normalisation.logSumNormalisation( scoreProfile, 0, 2*stop, normedProfile, 0 );
			
			System.arraycopy(normedProfile, 0, temp, 0, 2*stop );
			Arrays.sort( temp, 0, 2*stop );
			int e = 2*stop -1;
			double s = 0;
			while( e > 0 && (s+=temp[e]) < threshold ) {
				e--;
			}
			double t = temp[e];
			
			
			int k = 0;
			iList[0].clear();
			dList[0].clear();
			end.clear();
			tempPos.clear();
			for( e = 0; e < stop; e++ ) {
				if( normedProfile[e] >= t ) {
					normedProfile[k++] = scoreProfile[e];
					function[0].getLogScoreAndPartialDerivation( seq, start+e , iList[0], dList[0]); // motif gradient
					b.add( e );
					end.add( iList[0].length() );
				}
			}
			for( e = 0; e < stop; e++ ) {
				if( normedProfile[stop+e] >= t ) {
					normedProfile[k++] = scoreProfile[stop+e];
					try {
						function[0].getLogScoreAndPartialDerivation( seq.reverseComplement(), seq.getLength() - (start+e) - m , iList[0], dList[0]);
					} catch ( OperationNotSupportedException e1 ) {
						throw new RuntimeException( e1 );
					}
					b.add( -e-1 );
					end.add( iList[0].length() );
				}
			}
			//System.out.println("l: "+tempPos.length());
			
			if( add && start == 0 ) {
				toBeUsedHash.put( seq, b );
				//System.out.println(b+" "+seq);
			}
			
			componentScore[0] = logHiddenPotential[0] + Normalisation.logSumNormalisation( normedProfile, 0, k );
			
			for(k=0, counter = 0;k<b.length();k++){
				e = end.get( k );
				while( counter < e ){
					indices.add( iList[0].get( counter ) + paramRef[0] );
					partialDer.add( dList[0].get( counter++ ) * normedProfile[k] );
				}
			}
		} else {
			componentScore[0] = Double.NEGATIVE_INFINITY;
		}
		
		componentScore[function.length] = logHiddenPotential[function.length];

		double logScore = Normalisation.logSumNormalisation( componentScore, 0, componentScore.length, componentScore, 0 );
		
		partialDer.multiply( startIndex, partialDer.length(), componentScore[0] );
		
		// hiddenLambda
		for( int j = 0; j < 2; j++ )
		{
			indices.add( paramRef[1] + j );
			partialDer.add( componentScore[j] - (isNormalized()?hiddenPotential[j]:0) );
		}
		
		// bg for complete sequence
		return logScore + logP*(seq.getLength()-start);
	}

	@Override
	public String getInstanceName() {
		return getClass().getSimpleName() + "(" + function[0].getInstanceName() + ")";
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
			
			Sequence ts = sequence.getSubSequence( startpos, m.getLength() );
			
			double[] temp = new double[]{m.getLogScoreFor( ts ),
			                             m.getLogScoreFor( ts.reverseComplement() )};
			
			Normalisation.logSumNormalisation( temp );
			return temp;
		}
	}
	
	public void reset() {
		toBeUsedHash.clear();
	}

	public void setMotifModel( DifferentiableStatisticalModel bn ) {
		function[0] = bn;
		init( false );
	}
	
	public DifferentiableStatisticalModel getMotifModel(){
		return function[0];
	}
	
	public String toString(NumberFormat nf){
		StringBuffer buf = new StringBuffer();
		this.precomputeNorm();
		buf.append( "motif probability:"+ nf.format( Math.exp( partNorm[0] - norm ) )+"\n" );
		
		buf.append( function[0].toString() );
		
		return buf.toString();
	}
	
	public String toHtml(NumberFormat nf){
		StringBuffer buf = new StringBuffer();
		this.precomputeNorm();
		buf.append( "<p><strong>motif probability:</strong> "+ nf.format( Math.exp( partNorm[0] - norm ) )+"</p>" );
		
		if(function[0] instanceof BayesianNetworkDiffSM){
			buf.append( ((BayesianNetworkDiffSM)function[0]).toHtml(nf) );
		}else{
			buf.append( function[0].toString() );
		}
		
		return buf.toString();
	}
	
}
