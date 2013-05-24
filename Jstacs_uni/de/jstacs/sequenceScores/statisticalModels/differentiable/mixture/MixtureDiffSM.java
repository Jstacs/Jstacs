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

package de.jstacs.sequenceScores.statisticalModels.differentiable.mixture;

import java.text.NumberFormat;
import java.util.Arrays;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.motifDiscovery.MotifDiscoverer;
import de.jstacs.motifDiscovery.MutableMotifDiscoverer;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.VariableLengthDiffSM;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.DirichletMRGParams;

/**
 * This class implements a real mixture model.
 * 
 * @author Jens Keilwagen
 */
public class MixtureDiffSM extends AbstractMixtureDiffSM implements MutableMotifDiscoverer
{

	/**
	 * motifsRef[i] contains the number of motifs that are used in all components with index &lt; i.
	 */
	private int[] motifsRef;
	
	/**
	 * This constructor creates a new {@link MixtureDiffSM}. The first component determines the length of the sequences that can be modeled.
	 * 
	 * @param starts
	 *            the number of starts that should be done in an optimization
	 * @param plugIn
	 *            indicates whether the initial parameters for an optimization
	 *            should be related to the data or randomly drawn
	 * @param component
	 *            the {@link DifferentiableStatisticalModel}s
	 * 
	 * @throws CloneNotSupportedException
	 *             if an element of <code>component</code> could not be cloned
	 */
	public MixtureDiffSM( int starts, boolean plugIn, DifferentiableStatisticalModel... component ) throws CloneNotSupportedException {
		super( component[0].getLength(), starts, component.length, true, plugIn, component );
		for( int l, i = 0; i < component.length; i++ ) {
			l = component[i].getLength();
			if( l != 0 && length != l ) {
				throw new IllegalArgumentException( "The length of component " + i + " is " + l + " but should be " + length + "." );
			}
			if( !alphabets.checkConsistency( component[i].getAlphabetContainer() ) ) {
				throw new IllegalArgumentException( "The AlphabetContainer of component " + i + " is not suitable." );
			}
		}
		computeLogGammaSum();
		init();
	}

	private void init() {
		motifsRef = new int[function.length+1];
		for( int i = 0; i < function.length; i++ ) {
			motifsRef[i+1] = motifsRef[i];
			if( function[i] instanceof MotifDiscoverer ) {
				motifsRef[i+1] += ((MotifDiscoverer) function[i]).getNumberOfMotifs();
			}
		}
	}
	
	/**
	 * This is the constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link MixtureDiffSM} out of a
	 * {@link StringBuffer}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public MixtureDiffSM( StringBuffer xml ) throws NonParsableException {
		super( xml );
		init();
	}
	
	public MixtureDiffSM clone() throws CloneNotSupportedException {
		MixtureDiffSM clone = (MixtureDiffSM) super.clone();
		clone.init();
		return clone;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.AbstractMixtureDiffSM#getLogNormalizationConstantForComponent(int)
	 */
	@Override
	protected double getLogNormalizationConstantForComponent( int i ) {
		return function[i].getLogNormalizationConstant();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getLogPartialNormalizationConstant(int)
	 */
	@Override
	public double getLogPartialNormalizationConstant( int parameterIndex ) throws Exception {
		if( isNormalized() ) {
			return Double.NEGATIVE_INFINITY;
		} else {
			if( Double.isNaN( norm ) ) {
				precomputeNorm();
			}
			int[] ind = getIndices( parameterIndex );
			if( ind[0] == function.length ) {
				return partNorm[ind[1]];
			} else {
				return logHiddenPotential[ind[0]] + function[ind[0]].getLogPartialNormalizationConstant( ind[1] );
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.AbstractMixtureDiffSM#getHyperparameterForHiddenParameter(int)
	 */
	@Override
	public double getHyperparameterForHiddenParameter( int index ) {
		return function[index].getESS();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getEss()
	 */
	public double getESS() {
		double ess = 0;
		for( int i = 0; i < function.length; i++ ) {
			ess += function[i].getESS();
		}
		return ess;
	}

	
	private double[][] getRandomWeights( double[] originalWeights, int len ) {
		Arrays.fill( hiddenParameter, 0 );

		double[][] newWeights = new double[function.length][len];
		int i = 0, j = 0;
		double[] h = new double[this.getNumberOfComponents()];
		if( getESS() == 0 ) {
			Arrays.fill( h, 1 );
		} else {
			for( ; j < h.length; j++ ) {
				h[j] = getHyperparameterForHiddenParameter( j );
			}
		}
		DirichletMRGParams param = new DirichletMRGParams( h );
		double[] p = new double[h.length];
		double w = 1;
		while( i < newWeights[0].length ) {
			DirichletMRG.DEFAULT_INSTANCE.generate( p, 0, p.length, param );
			if( originalWeights != null ) {
				w = originalWeights[i];
			}
			for( j = 0; j < p.length; j++ ) {
				newWeights[j][i] = w * p[j];
				hiddenParameter[j] += newWeights[j][i];
			}
			i++;
		}
		computeHiddenParameter( hiddenParameter, true );
		return newWeights;
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.AbstractMixtureDiffSM#initializeUsingPlugIn(int, boolean, de.jstacs.data.DataSet[], double[][])
	 */
	@Override
	protected void initializeUsingPlugIn( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception {
		if( weights == null ) {
			weights = new double[data.length][];
		}
		double[] help = weights[index];
		double[][] newWeights = null;
		for( int r = 0; r < 3; r++ ) {
			if( r == 0 ) {
				newWeights = getRandomWeights( help, data[index].getNumberOfElements() );
			} else {
				for( int n = 0; n < data[index].getNumberOfElements(); n++ ) {
					fillComponentScores(data[index].getElementAt(n), 0);
					Normalisation.logSumNormalisation(componentScore);
					for( int k = 0; k < function.length; k++ ) {
						newWeights[k][n] = componentScore[k]*help[n];
					}
				}
			}			
			
			for( int i = 0; i < function.length; i++ ) {
				weights[index] = newWeights[i];
				function[i].initializeFunction( index, freeParams, data, weights );
			}
		}
		weights[index] = help;
	}
	
	/**
	 * Adjusts all hidden parameters including duration and mixture parameters according to the current values of the remaining parameters.
	 * 
	 * @param index the index of the class of this instance
	 * @param data the array of data for all classes
	 * @param weights the weights for all sequences in data
	 * @throws Exception thrown if the hidden parameters could not be adjusted
	 */
	public void adjustHiddenParameters( int index, DataSet[] data, double[][] weights ) throws Exception {
		if( weights == null ) {
			weights = new double[data.length][];
		}
		double[][] newWeights = getRandomWeights( weights[index], data[index].getNumberOfElements() );
		
		double[] help = weights[index];
		for( int i = 0; i < function.length; i++ ) {
			weights[index] = newWeights[i];
			if( function[i] instanceof MutableMotifDiscoverer ) {
				((MutableMotifDiscoverer) function[i]).adjustHiddenParameters( index, data, weights );
			}
		}
		weights[index] = help;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getInstanceName()
	 */
	public String getInstanceName() {
		String erg = "mixture(" + function[0].getInstanceName();
		for( int i = 1; i < function.length; i++ ) {
			erg += ", " + function[i].getInstanceName();
		}
		return erg + ")";
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.AbstractMixtureDiffSM#fillComponentScores(de.jstacs.data.Sequence, int)
	 */
	@Override
	protected void fillComponentScores( Sequence seq, int start ) {
		for( int i = 0; i < function.length; i++ ) {
			if( function[i] instanceof VariableLengthDiffSM ) {
				if( length != 0 ) {
					componentScore[i] = logHiddenPotential[i] + ((VariableLengthDiffSM)function[i]).getLogScoreFor( seq, start, start+length-1 );
				} else {
					componentScore[i] = logHiddenPotential[i] + ((VariableLengthDiffSM)function[i]).getLogScoreFor( seq, start, seq.getLength()-1 );
				}
			} else {
				componentScore[i] = logHiddenPotential[i] + function[i].getLogScoreFor( seq, start);
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLogScoreAndPartialDerivation
	 * (de.jstacs.data.Sequence, int, de.jstacs.utils.IntList,
	 * de.jstacs.utils.DoubleList)
	 */
	public double getLogScoreAndPartialDerivation( Sequence seq, int start, IntList indices, DoubleList partialDer ) {
		int i = 0, j = 0, k = paramRef.length - 1;
		k = paramRef[k] - paramRef[k - 1];
		for( ; i < function.length; i++ ) {
			iList[i].clear();
			dList[i].clear();
			if( function[i] instanceof VariableLengthDiffSM ) {
				if( length != 0 ) {
					componentScore[i] = logHiddenPotential[i] + ((VariableLengthDiffSM)function[i]).getLogScoreAndPartialDerivation( seq, start, start+length-1, iList[i], dList[i] );
				} else {
					componentScore[i] = logHiddenPotential[i] + ((VariableLengthDiffSM)function[i]).getLogScoreAndPartialDerivation( seq, start, seq.getLength()-1, iList[i], dList[i] );
				}
			} else {
				componentScore[i] = logHiddenPotential[i] + function[i].getLogScoreAndPartialDerivation( seq, start, iList[i], dList[i] );
			}
		}
		double logScore = Normalisation.logSumNormalisation( componentScore, 0, function.length, componentScore, 0 );
		for( i = 0; i < function.length; i++ ) {
			for( j = 0; j < iList[i].length(); j++ ) {
				indices.add( paramRef[i] + iList[i].get( j ) );
				partialDer.add( componentScore[i] * dList[i].get( j ) );
			}
		}
		for( j = 0; j < k; j++ ) {
			indices.add( paramRef[i] + j );
			partialDer.add( componentScore[j] - ( isNormalized() ? hiddenPotential[j] : 0 ) );
		}
		return logScore;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.SequenceScore#toString(java.text.NumberFormat)
	 */
	@Override
	public String toString( NumberFormat nf ) {
		if( Double.isNaN( norm ) ) {
			precomputeNorm();
		}
		
		StringBuffer erg = new StringBuffer( function.length * 1000 );
		for( int i = 0; i < function.length; i++ ) {
			erg.append( "p(" + i + ") = " + nf.format(isNormalized()?hiddenPotential[i]:Math.exp( partNorm[i] - norm )) + "\n"
						+ function[i].toString(nf)
						+ "\n" );
		}
		return erg.toString();
	}

	private int getComponentFor( int motif ) {
		int res = 0;
		while( motif >= motifsRef[res] ) {
			res++;
		}
		return res-1;
	}

	public void initializeMotif( int motifIndex, DataSet data, double[] weights ) throws Exception {
		int c = getComponentFor( motifIndex );
		if( function[c] instanceof MutableMotifDiscoverer ) {
			((MutableMotifDiscoverer) function[c]).initializeMotif( motifIndex - motifsRef[c], data, weights );
		} else {
			System.out.println( "WARNING: Not possible!" );
		}		
	}
	
	public void initializeMotifRandomly( int motif ) throws Exception {
		int c = getComponentFor( motif );
		if( function[c] instanceof MutableMotifDiscoverer ) {
			((MutableMotifDiscoverer) function[c]).initializeMotifRandomly( motif - motifsRef[c] );
		} else {
			System.out.println( "WARNING: Not possible!" );
		}
	}

	public boolean modifyMotif( int motifIndex, int offsetLeft, int offsetRight ) throws Exception {
		int c = getComponentFor( motifIndex );
		if( function[c] instanceof MutableMotifDiscoverer ) {
			boolean b = ((MutableMotifDiscoverer) function[c]).modifyMotif( motifIndex - motifsRef[c], offsetLeft, offsetRight );
			if( b ) {
				init( this.freeParams );
			}
			return b;
		} else {
			return false;
		}
	}

	public int getGlobalIndexOfMotifInComponent( int component, int motif ) {
		int res = motifsRef[component] + motif;
		if( res >= motifsRef[component+1] ) {
			throw new IndexOutOfBoundsException( "Component " + component + " has only " + (motifsRef[component+1]-motifsRef[component]) + " motifs." );
		} else {
			return res;
		}
	}

	public int getIndexOfMaximalComponentFor( Sequence sequence ) throws Exception {
		return this.getIndexOfMaximalComponentFor( sequence, 0 );
	}

	public int getMotifLength( int motif ) {
		int c = getComponentFor( motif );
		return ((MotifDiscoverer) function[c]).getMotifLength( motif - motifsRef[c] );
	}

	public int getNumberOfMotifs() {
		return motifsRef[function.length];
	}

	public int getNumberOfMotifsInComponent( int component ) {
		return motifsRef[component+1] - motifsRef[component];
	}

	public double[] getProfileOfScoresFor( int component, int motif, Sequence sequence, int startpos, KindOfProfile kind ) throws Exception {
		if( kind == KindOfProfile.UNNORMALIZED_JOINT ) {
			if( function[component] instanceof MotifDiscoverer ) {
				MotifDiscoverer md = (MotifDiscoverer) function[component];
				double[] prof = null, current;
				int c = md.getNumberOfComponents();
				for( int n, m, i = 0; i < c; i++ ) {
					n = md.getNumberOfMotifsInComponent( i );
					m = 0;
					while( m < n && md.getGlobalIndexOfMotifInComponent( i, m ) != motif ) {
						m++;
					}
					if( m < n ) {
						//component i contains the motif
						current = md.getProfileOfScoresFor( i, motif, sequence, startpos, kind );
						if( prof == null ) {
							prof = current;
						} else {
							for( int p = 0; p < prof.length; p++ ) {
								prof[p] = Normalisation.getLogSum( prof[p], current[p] );
							}
						}
					}
				}
				if( prof == null ) {
					throw new IllegalArgumentException();
				} else {
					return prof;
				}
			} else {
				throw new IllegalArgumentException();
			}
		} else {
			throw new OperationNotSupportedException( "Currently it is only allowed to used KindOfProfile.UNNORMALIZED_JOINT" );
		}
	}

	public double[] getStrandProbabilitiesFor( int component, int motif, Sequence sequence, int startpos ) throws Exception {
		//TODO
		return new double[]{0.5,0.5};
	}
}
