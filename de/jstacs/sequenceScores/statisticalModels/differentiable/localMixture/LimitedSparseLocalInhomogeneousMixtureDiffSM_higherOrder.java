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

package de.jstacs.sequenceScores.statisticalModels.differentiable.localMixture;

import java.text.NumberFormat;
import java.util.Arrays;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.motifDiscovery.Mutable;
import de.jstacs.sequenceScores.QuickScanningSequenceScore;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.random.DiMRGParams;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.DirichletMRGParams;
import de.jstacs.utils.random.FastDirichletMRGParams;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * Class for a sparse local inhomogeneous mixture (Slim) model.
 * The Slim model may be of higher order (i.e., considering multiple preceding positions) and
 * the number of preceding positions that are considered may be limited (LSlim model).
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder extends AbstractDifferentiableStatisticalModel implements Mutable, QuickScanningSequenceScore {

	private int order;
	private int distance;
	private double ess;
	private double[] localMixtureScore;
	private double[][] ancestorScore, depGrad;
	
	//parameters
	private double[][] componentMixtureParameters;
	private double[][][] ancestorMixtureParameters;
	private double[][][][] dependencyParameters;
	
	//logNorm constants
	private double[] componentMixtureLogNorm;
	private double[][] ancestorMixtureLogNorm;
	private double[][][] dependencyLogNorm;
	
	//potentials
	private double[][] componentMixturePotential;
	private double[][][] ancestorMixturePotential;
	private double[][][][] dependencyPotential;
	
	//index
	private int[] componentMixtureIndex;
	private int[][] ancestorMixtureIndex;
	private int[][][] dependencyIndex;
	private int numParameter;
	
	private double q;
	private double[] e;
	private double[][] logGamma;
	private PriorType type;
	private int[] seq;
	
	private double[][] scoreHash;
	
	/**
	 * The type of the prior used by the Slim model
	 * @author Jens Keilwagen
	 *
	 */
	public static enum PriorType {
		/**
		 * Bayesian Dirichlet likelihood equivalent, uniform prior
		 */
		BDeu,
		/**
		 * A simple mixture prior
		 */
		Simple_Mixture,
		/**
		 * A complex mixture prior
		 */
		Complex_Mixture
	}
	
	/**
	 * Creates a new Slim model with given number of components and maximum distance. Uses the BDeu prior.
	 * @param alphabets the alphabet of sequences the model is defined on
	 * @param length the length of the sequences that may be scores
	 * @param order the number of components, i.e., the number of preceding positions considered jointly
	 * @param distance the maximum distance of preceding positions considered
	 * @param ess the equivalent sample size
	 * @param q Parameter q of the mixture prior, ignored for BDeu prior
	 * @param t the type of the prior
	 * @throws IllegalArgumentException if the ess or other parameters are not allowed
	 */
	public LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder( AlphabetContainer alphabets, int length, int order, int distance, double ess, double q, PriorType t ) throws IllegalArgumentException {
		super(alphabets, length);
		if( !alphabets.isSimple() ) {
			throw new IllegalArgumentException("Check alphabet");
		}
		if( ess < 0 ) {
			throw new IllegalArgumentException("Check ess");
		}
		this.ess = ess;
		if( order < 0 ) {
			throw new IllegalArgumentException("Check number of local components");
		}
		this.order = order;
		this.distance = distance;
		
		this.type = t;
		if( q<=0 || q >= 1 ) {
			throw new IllegalArgumentException("Check q");
		}
		this.q = q;
		
		int a = (int) this.alphabets.getAlphabetLengthAt(0);
		dependencyParameters = new double[length][][][];
		componentMixtureParameters = new double[length][];
		ancestorMixtureParameters = new double[length][][];
		for( int l = 0; l < length; l++ ) {
			componentMixtureParameters[l] = new double[Math.min(l, order)+1];
			dependencyParameters[l] = new double[componentMixtureParameters[l].length][][];
			ancestorMixtureParameters[l] = new double[componentMixtureParameters[l].length][];
			int context=1, dist=1;
			for( int o = 0; o < componentMixtureParameters[l].length; o++ ) {
				dependencyParameters[l][o] = new double[context][a];
				if( o != 0 ) {
					dist=Math.min( l-o+1, distance );
				}
				ancestorMixtureParameters[l][o] = new double[dist];
				context*=a;
			}
		}

		initHash();
		init();
	}

	/**
	 * Creates a {@link LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder} model from its XML representation
	 * @param xml the XML representation
	 * @throws NonParsableException if XML could not be parsed
	 */
	public LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder(StringBuffer xml) throws NonParsableException {
		super(xml);
		initHash();
		init();
	}
	
	/**
	 * Returns the order of the Slim model
	 * @return the order
	 */
	public int getOrder() {
		return order;
	}

	/**
	 * Returns the maximum distance of preceding positions considered in the LSlim model.
	 * @return the maximum distance
	 */
	public int getDistance() {
		return distance;
	}
		
	private void init() {

		this.seq = new int[length];
		
		localMixtureScore = new double[order+1];
		ancestorScore = new double[order+1][length-1];
		e = new double[order+1];
		
		try {
			logGamma = ArrayHandler.clone( componentMixtureParameters );
			
			componentMixturePotential = ArrayHandler.clone( componentMixtureParameters );
			componentMixtureLogNorm = new double[length];
			componentMixtureIndex = new int[length];
			
			ancestorMixturePotential = ArrayHandler.clone( ancestorMixtureParameters );
			dependencyPotential = ArrayHandler.clone( dependencyParameters );
			
			dependencyLogNorm = new double[length][][];
			dependencyIndex = new int[length][][];
			ancestorMixtureLogNorm = new double[length][];
			ancestorMixtureIndex = new int[length][];
			for( int l = 0; l < length; l++ ) {
				dependencyLogNorm[l] = new double[componentMixtureParameters[l].length][];
				dependencyIndex[l] = new int[componentMixtureParameters[l].length][];
				for( int o = 0; o < componentMixtureParameters[l].length; o++ ) {
					dependencyLogNorm[l][o] = new double[dependencyParameters[l][o].length];
					dependencyIndex[l][o] = new int[dependencyParameters[l][o].length];
				}
				ancestorMixtureLogNorm[l] = new double[componentMixtureParameters[l].length];
				ancestorMixtureIndex[l] = new int[componentMixtureParameters[l].length];
			}

			int a = (int) alphabets.getAlphabetLengthAt(0);
			depGrad = new double[dependencyParameters[length-1][dependencyParameters[length-1].length-1].length][a];
		} catch( Exception e ) {
			//does not happen
			throw new RuntimeException(e);
		}
		
		precompute();
		
		numParameter = 0;
		for( int l = 0; l < length; l++ ) {
			componentMixtureIndex[l] = numParameter;
			numParameter += componentMixtureParameters[l].length;
			for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
				ancestorMixtureIndex[l][c] = numParameter;
				numParameter += ancestorMixtureParameters[l][c].length;
				for( int b = 0; b < dependencyParameters[l][c].length; b++ ) {
					dependencyIndex[l][c][b] = numParameter;
					numParameter += dependencyParameters[l][c][b].length;
				}
			}
		}
		
		if( type == PriorType.Complex_Mixture ) {
			precomputeLogGamma();
		}
	}
	
	private void initHash() {
		if(distance<8){
			scoreHash = new double[length][];
			for(int i=0;i<scoreHash.length;i++){
				int l = Math.min(distance+1, i+1);
				int dim = (int)Math.pow((int)alphabets.getAlphabetLengthAt(i),l);
				scoreHash[i] = new double[dim];
				Arrays.fill(scoreHash[i], Double.NaN);
			}
		}else{
			scoreHash = null;
		}
	}
	
	private void clearHash(){
		//System.out.println(this.hashCode()+" "+scoreHash);
		if(scoreHash != null){
			for(int i=0;i<scoreHash.length;i++){
				Arrays.fill(scoreHash[i], 1.0);
			}
		}
	}

	private void precompute() {
		for( int l = 0; l < length; l++ ) {
			componentMixtureLogNorm[l] = Normalisation.logSumNormalisation( componentMixtureParameters[l], 0, componentMixtureParameters[l].length, componentMixturePotential[l], 0 );
			for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
				ancestorMixtureLogNorm[l][c] = Normalisation.logSumNormalisation( ancestorMixtureParameters[l][c], 0, ancestorMixtureParameters[l][c].length, ancestorMixturePotential[l][c], 0 );
				for( int b = 0; b < dependencyParameters[l][c].length; b++ ) {
					dependencyLogNorm[l][c][b] = Normalisation.logSumNormalisation( dependencyParameters[l][c][b], 0, dependencyParameters[l][c][b].length, dependencyPotential[l][c][b], 0 );
				}
			}
		}
		clearHash();
	}

	public LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder clone() throws CloneNotSupportedException {
		LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder clone = (LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder) super.clone();
		clone.seq = seq.clone();

		
		clone.componentMixtureParameters = ArrayHandler.clone( componentMixtureParameters );
		clone.ancestorMixtureParameters = ArrayHandler.clone( ancestorMixtureParameters );
		clone.dependencyParameters = ArrayHandler.clone( dependencyParameters );
		
		clone.componentMixtureLogNorm = componentMixtureLogNorm.clone();
		clone.ancestorMixtureLogNorm = ArrayHandler.clone( ancestorMixtureLogNorm );
		clone.dependencyLogNorm = ArrayHandler.clone( dependencyLogNorm );
		
		clone.componentMixturePotential = ArrayHandler.clone( componentMixturePotential );
		clone.ancestorMixturePotential = ArrayHandler.clone( ancestorMixturePotential );
		clone.dependencyPotential = ArrayHandler.clone( dependencyPotential );
		
		clone.componentMixtureIndex = componentMixtureIndex.clone();
		clone.ancestorMixtureIndex = ArrayHandler.clone( ancestorMixtureIndex );
		clone.dependencyIndex = ArrayHandler.clone( dependencyIndex );
		
		clone.localMixtureScore = localMixtureScore.clone();
		clone.ancestorScore = ArrayHandler.clone( ancestorScore );
		clone.depGrad = ArrayHandler.clone( depGrad );
		clone.e = e.clone();
		clone.logGamma = ArrayHandler.clone( logGamma );

		if(scoreHash != null){
			//clone.scoreHash = ArrayHandler.clone(scoreHash);
		}
		
		clone.clearHash();
		return clone;
	}
	
	@Override
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getLogNormalizationConstant() {
		return 0;
	}

	@Override
	public double getLogPartialNormalizationConstant(int parameterIndex) throws Exception {
		return Double.NEGATIVE_INFINITY;
	}
	
	@Override
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		int a = (int) alphabets.getAlphabetLengthAt(0);
		double f = ess == 0 ? (a*(order+1)) : ess;//TODO
		
		for( int l = 0; l < length; l++ ) {
			//TODO
			if( type == PriorType.Complex_Mixture && componentMixtureParameters[l].length > 1 ) {
				int idx = (int) Math.floor( r.nextDouble() * componentMixtureParameters[l].length );
				double v = (1d-q)/(componentMixturePotential[l].length-1d);
				Arrays.fill( e, v*f );
				e[idx] = q*f;
			} else {
				Arrays.fill( e, f/componentMixtureParameters[l].length );
			}
			//draw( componentMixtureParameters[l], componentMixtureParameters[l].length, new DirichletMRGParams(e) );
			draw( componentMixtureParameters[l], componentMixtureParameters[l].length, new DirichletMRGParams(0,componentMixtureParameters[l].length, e) );
			for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
				int k=ancestorMixtureParameters[l][c].length;
				draw( ancestorMixtureParameters[l][c], k, new FastDirichletMRGParams( e[c]/k ) );
				FastDirichletMRGParams depPars = new FastDirichletMRGParams( e[c]/(dependencyParameters[l][c].length*a) );
				for( int b = 0; b < dependencyParameters[l][c].length; b++ ) {
					draw( dependencyParameters[l][c][b], a, depPars );	
				}
			}
		}
		precompute();
	}
	
	static void draw ( double[] array, int end, DiMRGParams pars ) {
		if( end > 1 ) {
			DirichletMRG.DEFAULT_INSTANCE.generateLog( array, 0, end, pars );
		}
	}
	
	@Override
	public void initializeFunction(int index, boolean freeParams, DataSet[] data, double[][] weights) throws Exception {
		
		//prepare
		for( int l = 0; l < length; l++ ) {
			Arrays.fill( componentMixtureParameters[l], Math.log( q/(componentMixtureParameters[l].length-1d) ) );//TODO
			componentMixtureParameters[l][0] = Math.log( 1d-q );//TODO
			for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
				//Arrays.fill( ancestorMixtureParameters[l][c], Math.exp( componentMixtureParameters[l][c] )/ancestorMixtureParameters[l][c].length );
				Arrays.fill( ancestorMixtureParameters[l][c], 0 );
				for( int b = 0; b < dependencyParameters[l][c].length; b++ ) {
					Arrays.fill( dependencyParameters[l][c][b], Math.exp( componentMixtureParameters[l][c] ) / dependencyParameters[l][c].length / dependencyParameters[l][c][b].length );
				}
			}
		}
		
		//count
		double w = 1;
		int a = alphabets.getAlphabetIndexForPosition(0), context;
		for( int i = 0; i < data[index].getNumberOfElements(); i++ ) {
			if( weights != null && weights[index] != null ) {
				w = weights[index][i];
			}
			Sequence s = data[index].getElementAt(i); 
			for( int l = 0; l < length; l++ ) {
				for( int c = 0; c < componentMixtureParameters[l].length; c++) {
					context = getOffset(s, l, c, dependencyParameters[0][0][0].length);
					for(int k=0;k<ancestorMixtureParameters[l][c].length;k++){
						context = next(context, s, 0, l, c, k );
						dependencyParameters[l][c][context][s.discreteVal( l )] += w;

					}
				}
			}
		}
		
		//divide and log
		for( int l = 0; l < length; l++ ) {
			for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
				for( int b = 0; b < dependencyParameters[l][c].length; b++ ) {
					Normalisation.sumNormalisation( dependencyParameters[l][c][b] );
					for( a = 0; a < dependencyParameters[l][c][b].length; a++){
						dependencyParameters[l][c][b][a] = Math.log( dependencyParameters[l][c][b][a]);
					}
				}
			}
		}
		
		precompute();
		//System.out.println(this);
	}
		
	private double unifGamma( double e, int anz ) {
		return anz * Gamma.logOfGamma( e/anz ) - Gamma.logOfGamma( e );
	}
	
	private double precomputePartLogGamma( int l ) {
		double lg = 0;
		for( int i = 0; i < componentMixtureParameters[l].length; i++ ) {
			lg += Gamma.logOfGamma( e[i] )*componentMixtureParameters[l][i];
		}
		lg -= Gamma.logOfGamma(this.ess);
		
		for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
			lg += unifGamma( e[c], ancestorMixtureParameters[l][c].length );
			for( int b = 0; b < dependencyParameters[l][c].length; b++ ) {
				lg += unifGamma( e[c]/dependencyParameters[l][c].length, dependencyParameters[l][c][b].length );
			}
		}
		return lg;
	}
	
	private void precomputeLogGamma() {
		//System.out.println("!!!!!!!! +" + q );
		double v, w;
		for( int l = 0; l < length; l++ ) {
			if( componentMixtureParameters[l].length == 1 ) {
				w = 1;
				v = 0;
			} else {
				w = q;
				v = (1d-q)/(componentMixturePotential[l].length-1d);
			}
			Arrays.fill(e, v*ess);
			for( int i = 0; i < componentMixtureParameters[l].length; i++ ) {
				e[i] = w*ess;
				logGamma[l][i] += precomputePartLogGamma(l);
				e[i] = v*ess;
			}
			//System.out.println( l + "\t" + Arrays.toString(logGamma[l]) );
		}
	}

	private static final double logPrior( double ess, double[] parameter, double logNorm ) {
		double res = 0, e = ess / parameter.length;
		for( int i = 0; i < parameter.length; i++ ) {
			res += e*parameter[i];
		}
		return res-ess*logNorm;
	}
	
	private void setLogPriorTerm( int l, int index, double logP ) {
		localMixtureScore[index]= logP;
		for( int i = 0; i < componentMixtureParameters[l].length; i++ ) {
			localMixtureScore[index] += e[i]*componentMixtureParameters[l][i];
		}
		localMixtureScore[index] -= this.ess*componentMixtureLogNorm[l];
		
		for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
			localMixtureScore[index] += logPrior( e[c], ancestorMixtureParameters[l][c], ancestorMixtureLogNorm[l][c] );
			for( int b = 0; b < dependencyParameters[l][c].length; b++ ) {
				localMixtureScore[index] += logPrior( e[c]/dependencyParameters[l][c].length, dependencyParameters[l][c][b], dependencyLogNorm[l][c][b] );	
			}
		}
	}
	
	@Override
	public double getLogPriorTerm() {
		double logPrior = 0, v, w, logP=0;
		for( int l = 0; l < length; l++ ) {
			if( componentMixtureParameters[l].length == 1 ) {
				w = 1;
				v = 0;
			} else {
				w = q;
				v = (1d-q)/(componentMixturePotential[l].length-1d);
			}
			switch( type ) {
				case BDeu:
					Arrays.fill(e, ess/componentMixtureParameters[l].length);
					setLogPriorTerm(l, 0, 1);
					logPrior += localMixtureScore[0];
					break;
				case Complex_Mixture:
					//logP =  Math.log(1d/componentMixtureParameters[l].length );
					Arrays.fill(e, v*ess);
					for( int i = 0; i < componentMixtureParameters[l].length; i++ ) {
						e[i] = w*ess;
						setLogPriorTerm(l, i, logP);
						localMixtureScore[i] += logGamma[l][i];
						//System.out.println( l + " " + i + "\t" + Arrays.toString(e) + "\t" + localMixtureScore[i] );
						e[i] = v*ess;
					}
					logPrior += Normalisation.getLogSum(0, componentMixtureParameters[l].length, localMixtureScore );
					break;
				case Simple_Mixture:
					for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
						localMixtureScore[c] = 0;
						for( int i = 0; i < componentMixtureParameters[l].length; i++ ) {
							localMixtureScore[c] += (i==c?w:v) * ess * componentMixtureParameters[l][i];	
						}
						localMixtureScore[c] -= ess * componentMixtureLogNorm[l];
					}
					logPrior += Normalisation.getLogSum( 0, componentMixtureParameters[l].length, localMixtureScore );
					double e  = ess / componentMixtureParameters[l].length;
					for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
						logPrior += logPrior( e, ancestorMixtureParameters[l][c], ancestorMixtureLogNorm[l][c] );
						for( int b = 0; b < dependencyParameters[l][c].length; b++ ) {
							logPrior += logPrior( e/dependencyParameters[l][c].length, dependencyParameters[l][c][b], dependencyLogNorm[l][c][b] );	
						}
					}
					break;
			}
		}
		return logPrior;
	}
	
	private static final void addGrad( double ess, double[] potential, double[] grad, int start, double w ) {
		double e = ess / potential.length;
		for( int i = 0; i < potential.length; i++ ) {
			grad[start+i] += w*(e - ess*potential[i]);
		}
	}

	private void addGradientPartOfLogPriorTerm(double[] grad, int start, int l, double w ) throws Exception {
		for( int i = 0; i < componentMixtureParameters[l].length; i++ ) {
			grad[start+componentMixtureIndex[l]+i] += w*(e[i] - this.ess * componentMixturePotential[l][i]);
		}

		for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
			addGrad( e[c], ancestorMixturePotential[l][c], grad, start+ancestorMixtureIndex[l][c], w );
			for( int b = 0; b < dependencyParameters[l][c].length; b++ ) {
				addGrad( e[c]/dependencyParameters[l][c].length, dependencyPotential[l][c][b], grad, start+dependencyIndex[l][c][b], w );
			}
		}
	}
	
	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int start) throws Exception {
		double v, w, logP = 0;
		for( int l = 0; l < length; l++ ) {
			if( componentMixtureParameters[l].length == 1 ) {
				w = 1;
				v = 0;
			} else {
				w = q;
				v = (1d-q)/(componentMixturePotential[l].length-1d);
			}
			switch( type ) {
				case BDeu:
					Arrays.fill(e, ess/componentMixtureParameters[l].length);
					addGradientPartOfLogPriorTerm(grad, start, l, 1);
					break;
				case Complex_Mixture:
					//logP =  Math.log(1d/componentMixtureParameters[l].length );
					Arrays.fill(e, v*ess);
					for( int i = 0; i < componentMixtureParameters[l].length; i++ ) {
						e[i] = w*ess;
						setLogPriorTerm(l, i, logP);
						e[i] = v*ess;
					}
					Normalisation.logSumNormalisation( localMixtureScore, 0, componentMixtureParameters[l].length );
					for( int i = 0; i < componentMixtureParameters[l].length; i++ ) {
						e[i] = w*ess;
						addGradientPartOfLogPriorTerm(grad, start, l, localMixtureScore[i]);
						e[i] = v*ess;
					}
					break;
				case Simple_Mixture:
					for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
						localMixtureScore[c] = 0;
						for( int i = 0; i < componentMixtureParameters[l].length; i++ ) {
							localMixtureScore[c] += (i==c?w:v) * ess * componentMixtureParameters[l][i];	
						}
						localMixtureScore[c] -= ess * componentMixtureLogNorm[l];
					}
					Normalisation.logSumNormalisation( localMixtureScore, 0, componentMixtureParameters[l].length );
					for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
						for( int i = 0; i < componentMixtureParameters[l].length; i++ ) {
							grad[start+componentMixtureIndex[l]+c] += localMixtureScore[i] * (i==c?w:v) * ess;	
						}
						grad[start+componentMixtureIndex[l]+c] -= ess*componentMixturePotential[l][c]; 
					}
					double e  = ess / componentMixtureParameters[l].length;
					for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
						addGrad( e, ancestorMixturePotential[l][c], grad, start+ancestorMixtureIndex[l][c], 1 );
						for( int b = 0; b < dependencyParameters[l][c].length; b++ ) {
							addGrad( e/dependencyParameters[l][c].length, dependencyPotential[l][c][b], grad, start+dependencyIndex[l][c][b], 1 );
						}
					}
					break;
			}
		}
	}

	@Override
	public double getESS() {
		return ess;
	}
	
	
	
	private int getOffset( Sequence seq, int start, int order, int a ) {
		int offset = 0;
		for( int o = 1; o < order; o++ ) {
			offset = offset*a + seq.discreteVal(start-o);
		}
		return offset;
	}
	
	private int getOffset( int[] seq, int start, int order, int a ) {
		int offset = 0;
		for( int o = 1; o < order; o++ ) {
			offset = offset*a +seq[start-o];
		}
		return offset;
	}
	
	private int next( int context, Sequence seq, int start, int l, int c, int m ) {
		return (context*dependencyParameters[l][c][0].length + seq.discreteVal(start+l-c-m))%dependencyParameters[l][c].length;
	}
	
	private int next( int context, int[] seq, int start, int l, int c, int m ) {
		return (context*dependencyParameters[l][c][0].length + seq[start+l-c-m])%dependencyParameters[l][c].length;
	}
	
	public boolean[][] getInfixFilter( int kmer, double thresh, int... start ) {
		if( order > 1 ) {
			throw new RuntimeException("Not implemented");
		}
		
		double[][] val = getCum_Complex(kmer);
		double[] cum = val[0];
		double[] revCum = val[1];
		
		int a = dependencyParameters[0][0][0].length;
		boolean[][] use = new boolean[start.length][(int)Math.pow(a, kmer)];

		//more tricky (faster)
		int maxSymbol = a-1;
		int z=0;
		int[] seq = new int[kmer + 1];
		double[][] prefixScore = new double[start.length][kmer+1];
		int l=0;
		double[] best=new double[start.length];
		Arrays.fill( best, Double.NEGATIVE_INFINITY );

		do {
			for( int s = 0; s  < start.length; s++ ) {
				getInfixScores(s, start[s], l, kmer, seq, prefixScore, null);
				if( prefixScore[s][kmer]> best[s] ) {
					best[s] = prefixScore[s][kmer];
				}
				use[s][z]= (cum[start[s]] + prefixScore[s][kmer]+revCum[start[s]+kmer]) >= thresh;
			}
			
			z++;
			l = kmer-1-next(seq, maxSymbol);
		} while( seq[kmer]==0 );
/*		
System.out.println("cum "+ Arrays.toString(cum) );
System.out.println("revCum "+ Arrays.toString(revCum) );
*/
		return use;
	}
	
	protected double[][] getCum_Naive( int kmer ) {
		int len = getLength();
		double[] max = new double[len];
		double[] cum = new double[len+1];
		double[] revCum = new double[len+1];
		double last =0;
		for( int l = 0; l < max.length; l++ ) {
			//compute maximum for position l independent of all other positions
			max[l] = Double.NEGATIVE_INFINITY;
			for( int a = 0; a < dependencyPotential[l][0][0].length; a++ ) {
				double v=0;
				for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
					int an = 0;
					for( int i = 1; i < dependencyParameters[l][c].length; i++ ) {
						if( dependencyPotential[l][c][an][a] < dependencyPotential[l][c][i][a] ) {
							an=i;
						}
					}
					localMixtureScore[c] = componentMixtureParameters[l][c] - componentMixtureLogNorm[l] + dependencyParameters[l][c][an][a] - dependencyLogNorm[l][c][an];
				}
				v=Normalisation.getLogSum(0, componentMixtureParameters[l].length, localMixtureScore );
				if( v > max[l] ) {
					max[l] = v;
				}
			}
			
			cum[l+1] = last + max[l];
			last = cum[l+1];
		}
		last = 0;
		for( int l = len-1; l>=0; l-- ) {
			revCum[l] = last+max[l];
			last = revCum[l];
		}
		return new double[][]{cum, revCum};
	}
	
	protected double[][] getCum_Complex( int kmer ) {
		int start = getLength()-kmer+1;
		int len = getLength();
		double[][] max = new double[kmer][]; // max[k][l] = maximum score of a infix of length (k+1) at position l
		for( int k = 0; k < max.length; k++ ) {
			max[k] = new double[len-k];
			Arrays.fill( max[k],  Double.NEGATIVE_INFINITY );
		}
		
		int a = dependencyParameters[0][0][0].length, i;

		//fill max table
		double[][] prefixScore = new double[len][kmer+1];
		int maxSymbol = a-1;
		int[] seq = new int[kmer + 1];
		
		// compute all max values where a full kmer if possible
		int l=0;
		do {
			for( int s = 0; s < start; s++ ) {
				getInfixScores(s, s, l, kmer, seq, prefixScore, max);
			}			
			// increment the sequence seq and save the position of the modification (l)
			l = kmer-1-next(seq, maxSymbol);
		} while( seq[kmer]==0 );

		//fill missing values of max array (g < k)
		for( int s = start; s < len; s++ ) {
			int gmer = len-s;
			Arrays.fill(seq, 0);
			l=0;
			do {
				getInfixScores(s, s, l, gmer, seq, prefixScore, max);
				l = gmer-1-next(seq, maxSymbol);
			} while( seq[gmer]==0 );
		}
		
		
		//DP algorithms for cumulatives => we like to find the tightest bound
		
		//cumulative
		double[] cum = new double[len+1];
		Arrays.fill(cum, Double.POSITIVE_INFINITY);
		cum[0] = 0;
		for( i = 0; i < len; i++ ) {
			for( int j = Math.min(i, kmer-1); j>=0; j-- ) {
				if( cum[i+1] > cum[i-j]+max[j][i-j] ) {
					cum[i+1] = cum[i-j]+max[j][i-j];
				}
			}
		}
		//rev. cumulative
		double[] revCum = new double[len+1];
		Arrays.fill(revCum, Double.POSITIVE_INFINITY);
		revCum[len] = 0;
		for( i = len-1; i >= 0; i-- ) {
			for( int j = Math.min(len-1-i, kmer-1); j>=0; j-- ) {
				if( revCum[i] > revCum[i+1+j]+max[j][i] ) {
					revCum[i] = revCum[i+1+j]+max[j][i];
				}
			}
		}
		return new double[][]{cum, revCum};
	}

	protected void getInfixScores(int s, int start, int l, int kmer, int[] seq, double[][] prefixScore, double[][] max) {
		while( l < kmer ) {
			int pos = l+start; // position in (L)slim model
			int ll = kmer-1-l; //position in the sequence seq
			int current = seq[ll];
			
			//unconditional score
			localMixtureScore[0] = componentMixtureParameters[pos][0] - componentMixtureLogNorm[pos] + dependencyParameters[pos][0][0][current] - dependencyLogNorm[pos][0][0];

			//conditional scores
			for( int c = 1; c < componentMixtureParameters[pos].length; c++ ) {
				int g=ancestorMixtureParameters[pos][c].length;
				for( int m = 0; m < g; m++ ) {
					int an;
					if( ll+c+m < kmer ) {
						//context known
						an = seq[ll+c+m];
					} else {
						//context not known
						an=0;
						for( int i = 1; i < dependencyParameters[pos][c].length; i++ ) {
							if( dependencyPotential[pos][c][an][current] < dependencyPotential[pos][c][i][current] ) {
								an=i;
							}
						}
					}
					ancestorScore[c][m] = ancestorMixtureParameters[pos][c][m] - ancestorMixtureLogNorm[pos][c] + dependencyParameters[pos][c][an][current] - dependencyLogNorm[pos][c][an];
				}
				localMixtureScore[c] = componentMixtureParameters[pos][c] - componentMixtureLogNorm[pos] + Normalisation.getLogSum(0, g, ancestorScore[c]);
			}
			
			//complete score
			double score = Normalisation.getLogSum(0, componentMixtureParameters[pos].length, localMixtureScore );
			
			prefixScore[s][l+1]=prefixScore[s][l]+score;
			if( max!= null && prefixScore[s][l+1] > max[l][s] ) {
				max[l][s] = prefixScore[s][l+1];
			}				
			
			l++;
		}
	}
	
	private static int next( int[] seq, int maxSymbol ) {
		int i=0;
		while( seq[i] == maxSymbol ) {
			seq[i++] = 0;
		}
		seq[i]++;
		return i;
	}	
	
	@Override
	public double getLogScoreFor(Sequence seq2, int start) {
		for(int i=0;i<this.seq.length;i++){
			this.seq[i] = seq2.discreteVal(start+i);
		}
		int al = (int)alphabets.getAlphabetLengthAt(0);
		int pw = (int)Math.pow(al, distance);
		double score = 0;
		int idx = 0;
		for( int l = 0; l < length; l++ ) {
			if(l>distance){
				idx %= pw;
			}
			idx = idx*al + seq[l];
			
			double hscore = 1;
			if(scoreHash != null){
				hscore = scoreHash[l][idx];
			}
			/*if(hscore > 0 && hscore != 1.0){
				System.out.println(hscore);
			}*/
			if(hscore > 0){
				int current = seq[l];
				localMixtureScore[0] = componentMixtureParameters[l][0] - componentMixtureLogNorm[l] + dependencyParameters[l][0][0][current] - dependencyLogNorm[l][0][0];
				for( int c = 1; c < componentMixtureParameters[l].length; c++ ) {
					int k=ancestorMixtureParameters[l][c].length;
					int an=getOffset(seq, l, c, dependencyParameters[0][0][0].length);

					int l1 = dependencyParameters[l][c][0].length;
					int l2 = dependencyParameters[l][c].length;

					for( int m = 0; m < k; m++ ) {
						//an = next(an, seq, start, l, c, m );
						an = (an*l1 + seq[l-c-m])%l2;
						ancestorScore[c][m] = ancestorMixtureParameters[l][c][m] - ancestorMixtureLogNorm[l][c] + dependencyParameters[l][c][an][current] - dependencyLogNorm[l][c][an]; 
					}
					localMixtureScore[c] = componentMixtureParameters[l][c] - componentMixtureLogNorm[l] + Normalisation.getLogSum(0, k, ancestorScore[c]);
				}
				hscore = Normalisation.getLogSum(0, componentMixtureParameters[l].length, localMixtureScore );
				if(scoreHash != null){
					scoreHash[l][idx] = hscore;
				}
				
			}
			score += hscore;
			
		}

		return score;
	}

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList partialDer) {
		double score = 0;
		for( int l = 0; l < length; l++ ) {
			int current = seq.discreteVal(start+l);
			localMixtureScore[0] = componentMixtureParameters[l][0] - componentMixtureLogNorm[l] + dependencyParameters[l][0][0][current] - dependencyLogNorm[l][0][0];
			for( int c = 1; c < componentMixtureParameters[l].length; c++ ) {
				int k=ancestorMixtureParameters[l][c].length;
				int an=getOffset(seq, start+l, c, dependencyParameters[0][0][0].length);
				for( int m = 0; m < k; m++ ) {
					an = next(an, seq, start, l, c, m );
					ancestorScore[c][m] = ancestorMixtureParameters[l][c][m] - ancestorMixtureLogNorm[l][c] + dependencyParameters[l][c][an][current] - dependencyLogNorm[l][c][an]; 
				}
				localMixtureScore[c] = componentMixtureParameters[l][c] - componentMixtureLogNorm[l] + Normalisation.logSumNormalisation(ancestorScore[c], 0, k);
			}
			
			
			//gradient
			score += Normalisation.logSumNormalisation( localMixtureScore, 0, componentMixtureParameters[l].length );
			for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
				indices.add( componentMixtureIndex[l]+c );
				partialDer.add( localMixtureScore[c] - componentMixturePotential[l][c] );
			}
			for( int a = 0; a < dependencyParameters[l][0][0].length; a++ ) {
				indices.add( dependencyIndex[l][0][0]+a );
				partialDer.add( localMixtureScore[0]*( (a==current?1:0) - dependencyPotential[l][0][0][a] ) );
			}
			for( int c = 1; c < componentMixtureParameters[l].length; c++ ) {
				
				for( int a = 0; a < dependencyParameters[l][c].length; a++ ) {
					Arrays.fill( depGrad[a], 0 );
				}
				int k=ancestorMixtureParameters[l][c].length;
				
				int an=getOffset(seq, start+l, c, dependencyParameters[0][0][0].length);
				for( int m = 0; m < k; m++ ) {
					indices.add( ancestorMixtureIndex[l][c]+m );
					partialDer.add( localMixtureScore[c] * ( ancestorScore[c][m] - ancestorMixturePotential[l][c][m] ) );
					
					an = next(an, seq, start, l, c, m );
					for( int a = 0; a < depGrad[an].length; a++ ) {
						depGrad[an][a] += ancestorScore[c][m] * ((a==current?1:0)-dependencyPotential[l][c][an][a]);
					}
				}
				
				for( int a = 0; a < dependencyParameters[l][c].length; a++ ) {
					for( int b = 0; b < dependencyParameters[l][c][a].length; b++ ) {
						if( depGrad[a][b] != 0 ) {
							indices.add( dependencyIndex[l][c][a]+b );
							partialDer.add( localMixtureScore[c] * depGrad[a][b] );
						}
					}
				}
			}
		}
		return score;
	}
	
	//TODO remove
	/*public void checkGrad( IntList indices, DoubleList partDer ) {
		double[] grad = new double[numParameter];
		for( int i = 0; i < indices.length(); i++ ) {
			grad[indices.get(i)] += partDer.get(i);
		}
		System.out.println( Arrays.toString(grad) );
		for( int l = 0; l < length; l++ ) {
			System.out.println( "comMix " + l + " " + ToolBox.sum( componentMixtureIndex[l], componentMixtureIndex[l]+componentMixtureParameters[l].length, grad) );
			for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
				System.out.println( "anMix " + l + " " + c + " " + ToolBox.sum( ancestorMixtureIndex[l][c], ancestorMixtureIndex[l][c]+ancestorMixtureParameters[l][c].length, grad) );
				for( int b = 0; b < dependencyParameters[l][c].length; b++ ) {
					System.out.println( "dep " + l + " " + c + " " + b + " " + ToolBox.sum( dependencyIndex[l][c][b], dependencyIndex[l][c][b]+dependencyParameters[l][c][b].length, grad) );
				}
			}
		}
	}*/

	@Override
	public int getNumberOfParameters() {
		//System.out.println("numPars: "+numParameter);
		return numParameter;
	}

	private void get( int start, double[] parameter, double[] array ) {
		for( int i = 0; i < parameter.length; i++ ) {
			array[start+i] = parameter[i];
		}
	}
	
	@Override
	public double[] getCurrentParameterValues() throws Exception {
		double[] res = new double[numParameter];
		//System.out.println("returned: "+res.length);
		for( int l = 0; l < length; l++ ) {
			get( componentMixtureIndex[l], componentMixtureParameters[l], res );
			for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
				get( ancestorMixtureIndex[l][c], ancestorMixtureParameters[l][c], res );
				for( int b = 0; b < dependencyParameters[l][c].length; b++ ) {
					get( dependencyIndex[l][c][b], dependencyParameters[l][c][b], res );
				}
			}
		}
		return res;
	}

	private void set( int start, double[] parameter, double[] array ) {
		for( int i = 0; i < parameter.length; i++ ) {
			parameter[i] = array[start+i];
		}
	}
	
	/**
	 * Sets the (conditional) probability parameters at a specific position and sets the mixture parameters
	 * (largely) to the unconditional PWM component.
	 * @param position the position
	 * @param pars the new parameters at this position
	 */
	public void set(int position, double[] pars){
		for(int i=0;i<dependencyParameters[position].length;i++){
			for(int j=0;j<dependencyParameters[position][i].length;j++){
				System.arraycopy( pars, 0, dependencyParameters[position][i][j], 0, pars.length );
			}
		}
		Arrays.fill( componentMixtureParameters[position], -1000);
		componentMixtureParameters[position][0] = 0;
		
		precompute();
	}
	
	@Override
	public void setParameters(double[] params, int start) {
		for( int l = 0; l < length; l++ ) {
			set( start+componentMixtureIndex[l], componentMixtureParameters[l], params);
			for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
				set( start+ancestorMixtureIndex[l][c], ancestorMixtureParameters[l][c], params);
				for( int b = 0; b < dependencyParameters[l][c].length; b++ ) {
					set( start+dependencyIndex[l][c][b], dependencyParameters[l][c][b], params);
				}
			}
		}
		precompute();
	}

	@Override
	public String getInstanceName() {
		return "slim(" + order + ", " + distance + ", " + type + " " + ess + " " + q + ")";
	}

	@Override
	public boolean isInitialized() {
		return true;
	}

	private static final String XML_TAG = "SLIM";
	
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, alphabets, "alphabets");
		XMLParser.appendObjectWithTags(xml, length, "length");
		XMLParser.appendObjectWithTags(xml, order, "order");
		XMLParser.appendObjectWithTags( xml, distance, "distance" );
		XMLParser.appendObjectWithTags(xml, ess, "ess");
		
		XMLParser.appendObjectWithTags(xml, componentMixtureParameters, "componentMixtureParameters");
		XMLParser.appendObjectWithTags(xml, ancestorMixtureParameters, "ancestorMixtureParameters");
		XMLParser.appendObjectWithTags(xml, dependencyParameters, "dependencyParameters");
		
		XMLParser.appendObjectWithTags(xml, type, "priorType");
		XMLParser.appendObjectWithTags(xml, q, "q");
		
		XMLParser.addTags(xml, XML_TAG);
		return xml;
	}

	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		xml = XMLParser.extractForTag(xml, XML_TAG);
		
		alphabets = (AlphabetContainer) XMLParser.extractObjectForTags(xml, "alphabets");
		length = (Integer) XMLParser.extractObjectForTags(xml, "length");
		order = (Integer) XMLParser.extractObjectForTags(xml, "order");
		distance = (Integer) XMLParser.extractObjectForTags( xml, "distance" );
		ess = (Double) XMLParser.extractObjectForTags(xml, "ess");
		
		componentMixtureParameters = (double[][]) XMLParser.extractObjectForTags(xml, "componentMixtureParameters");
		ancestorMixtureParameters = (double[][][]) XMLParser.extractObjectForTags(xml, "ancestorMixtureParameters");
		dependencyParameters = (double[][][][]) XMLParser.extractObjectForTags(xml, "dependencyParameters");
		
		type = (PriorType) XMLParser.extractObjectForTags(xml, "priorType");
		q = (Double) XMLParser.extractObjectForTags(xml, "q");
	}
	
	public String toString( NumberFormat nf ) {
		StringBuffer res = new StringBuffer();
		DiscreteAlphabet abc = (DiscreteAlphabet) alphabets.getAlphabetAt(0);
		for( int l = 0; l < length; l++ ) {
			//TODO
			for( int c = 0; c < componentMixtureParameters[l].length; c++ ) {
				res.append( "P_" + l + "(c="+c+")=" + nf.format(componentMixturePotential[l][c]) + "\n" );
				if( c!= 0 ) {
					res.append( "P_" + l + "(p|c="+c+")=" + ToolBox.toString(ancestorMixturePotential[l][c],nf) + "\n" );
				}
				for( int b = 0; b < dependencyParameters[l][c].length; b++ ) {
					res.append( "P_" + l + "(x|" );
					if( c > 0 ) {
						int h = b;
						res.append( "y=" );
						for( int j =0; j < c; j++ ) {
							res.append( abc.getSymbolAt(h%dependencyParameters[0][0][0].length) );
							h /= dependencyParameters[0][0][0].length;
						}
						res.append( ",");
					}
					res.append("c="+c+")=" + ToolBox.toString(dependencyPotential[l][c][b],nf) + "\n" );
				}
				res.append("\n");
			}
			res.append("\n");
		}
		return res.toString() + getGraphviz();
	}
	
	/**
	 * Returns the conditional probabilities for the specified component.
	 * @param component the component
	 * @return the conditional probabilities
	 * @throws CloneNotSupportedException if the internal probabilities could not be cloned
	 */
	public double[][][] getConditionalProbabilities(int component) throws CloneNotSupportedException{
		double[][][] condProbs = new double[dependencyPotential.length][][];
		for(int i=0;i<condProbs.length;i++){
			if(dependencyPotential[i].length > component){
				condProbs[i] = ArrayHandler.clone( dependencyPotential[i][component] );
			}else{
				condProbs[i] = new double[0][0];
			}
		}
		return condProbs;
	}
	
	/**
	 * Returns the unconditional, normalized (PWM) probabilities of this Slim model
	 * @return the PWM probabilities
	 * @throws CloneNotSupportedException if the internal parameters could not be cloned
	 */
	public double[][] getPWMParameters() throws CloneNotSupportedException{
		double[][] pwm = new double[length][];
		for( int l = 0; l < length; l++ ) {
			pwm[l] = dependencyParameters[l][0][0].clone();
			Normalisation.logSumNormalisation( pwm[l] );
		}
		return pwm;
	}
	
	/**
	 * Returns the probabilities of the mixture components.
	 * @return the probabilities of the mixture components
	 * @throws CloneNotSupportedException if the internal parameters could not be cloned
	 */
	public double[][] getMixtureProbabilities() throws CloneNotSupportedException{
		return ArrayHandler.clone( componentMixturePotential );
	}
	
	/**
	 * Returns the probabilities that the preceding positions considered are used as context.
	 * @param component the component considered
	 * @return the probabilities of the preceding positions
	 */
	public double[][] getAncestorProbabilities(int component){
		double[][] ancProb = new double[ancestorMixturePotential.length][];
		for(int i=0;i<ancProb.length;i++){
			if(ancestorMixturePotential[i].length > component){
				ancProb[i] = ancestorMixturePotential[i][component].clone();
			}else{
				ancProb[i] = new double[0];
			}
		}
		return ancProb;
	}
	
	/*
	public String getGraphviz2( double s ) {
		StringBuffer res = new StringBuffer();
		for( int l = 0; l < length; l++ ) {
			for( int c = 1; c < componentMixtureParameters[l].length; c++ ) {
				for( int t=0; t < ancestorMixturePotential[l][c-1].length; t++ ) {
					double r = -100;
					for( int y=0; r<s && y < dependencyPotential[l][c-1].length; y++ ) {
						for( int x=0; r<s && x < dependencyPotential[l][c-1][y].length; x++ ) {
							//System.out.println("!"+l+"\t"+c+"\t"+t+"\t"+y+"\t"+x);
							r = componentMixturePotential[l][c] * ancestorMixturePotential[l][c-1][t]*dependencyPotential[l][c-1][y][x] 
									- componentMixturePotential[l][0] * symbolPotential[l][x];
						}
					}
					System.out.println(l + "\t"+t +"\t"+ r);
					if( r >= s ) {
						res.append( (l-t-1) + " -> " + l + "[label=" + r+ "]\n" );//TODO weight
					}
				}
			}
		}
		return res.toString();
	}*/
	
	/**
	 * Returns a Graphviz (dot) representation of the Slim model.
	 * @return the Graphviz representation
	 */
	public String getGraphviz() {
		StringBuffer res = new StringBuffer();
		for( int l = 0; l < length; l++ ) {
			for( int c = 1; c < componentMixtureParameters[l].length; c++ ) {
				for( int p = 0; p < ancestorMixturePotential[l][c].length; p++ ) {
					String w = Integer.toHexString( (int) (256*(1-componentMixturePotential[l][c]*ancestorMixturePotential[l][c][p])) );
					if( w.length() == 1 ) {
						w = "0" + w;
					}
					res.append( (l-p-1) + " -> " + l + "["
							//+"weight=" + nf.format(componentMixturePotential[l][c]*ancestorMixturePotential[l][c][p])+ ", "
							+ "color=\"#"+w+""+w+""+w+"\"]\n" );
				}
			}
		}
		return res.toString();
	}
	
	public int getNumberOfRecommendedStarts() {
		return 1;//TODO
	}

	@Override
	public boolean modify(int offsetLeft, int offsetRight) {
		if( offsetLeft != 0 || offsetRight != 0 ) {
			int oldLength = length;
			length = oldLength - offsetLeft + offsetRight;
			boolean modLength = oldLength != length; 
			int a = (int) alphabets.getAlphabetLengthAt(0);
			double[][][][] newDependencyParameters = new double[length][][][];
			double[][] newComponentMixtureParameters = new double[length][];
			double[][][] newAncestorMixtureParameters = new double[length][][];
						
			for( int oldIndex = offsetLeft, l = 0; l < length; l++, oldIndex++ ) {
				int m = Math.min(l, order)+1;
				newAncestorMixtureParameters[l] = new double[m][];
				
				if( oldIndex >= 0 && oldIndex < oldLength ) {
					int n = Math.min(m, componentMixtureParameters[oldIndex].length ); //n>0
					
					//mixtures and emissions
					if( m == componentMixtureParameters[oldIndex].length ) {
						//same order => copy completely
						newDependencyParameters[l] = dependencyParameters[oldIndex];
						newComponentMixtureParameters[l] = componentMixtureParameters[oldIndex];
					} else {
						//different order
						newDependencyParameters[l] = new double[m][][];
						newComponentMixtureParameters[l] = new double[m];
						//copy common part
						for( int i = 0; i < n; i++ ) {
							double sum = ToolBox.sum(0, n, componentMixturePotential[oldIndex]), f = n/(double)(m);
							newComponentMixtureParameters[l][i] = Math.log( componentMixturePotential[oldIndex][i]/sum * f );
							newDependencyParameters[l][i] = dependencyParameters[oldIndex][i];
						}
						//fill with last
						for( int i = n; i < m; i++ ) {
							//newDependencyParameters[l][i] = new double[newDependencyParameters[l][i-1].length*4][a];//TODO has been 4
							newDependencyParameters[l][i] = new double[newDependencyParameters[l][i-1].length*a][a];
							for( int b = 0, j = 0; j < newDependencyParameters[l][i-1].length; j++ ) {
								for( int c = 0; c < a; c++, b++ ) {
									System.arraycopy(newDependencyParameters[l][i-1][j], 0, newDependencyParameters[l][i][b], 0, newDependencyParameters[l][i-1][j].length);
								}
							}
						}
						Arrays.fill( newComponentMixtureParameters[l], n, m, -Math.log(m) );
					}

					//ancestorMixtureParameters
					int dist;
					newAncestorMixtureParameters[l][0] = new double[1];
					for( int o = 1; o < newAncestorMixtureParameters[l].length; o++ ) {
						dist = Math.min( l-o+1, distance );
						newAncestorMixtureParameters[l][o] = new double[dist];
						Arrays.fill( newAncestorMixtureParameters[l][o], -Math.log(dist) );
						if( o < ancestorMixtureParameters[oldIndex].length ) {
							int h = Math.min( dist, ancestorMixtureParameters[oldIndex][o].length );
							double sum = ToolBox.sum(0, h, ancestorMixturePotential[oldIndex][o]), f = h/(double)dist;
							for( int j = 0; j < h; j++ ) {
								newAncestorMixtureParameters[l][o][j] = Math.log( ancestorMixturePotential[oldIndex][o][j]/sum * f );
							}
						}
					}
				} else {
					//uniform
					newComponentMixtureParameters[l] = new double[m];
					newDependencyParameters[l] = new double[newComponentMixtureParameters[l].length][][];
					newAncestorMixtureParameters[l] = new double[newComponentMixtureParameters[l].length][];
					int context=1, dist = 1;
					for( int o = 0; o < newComponentMixtureParameters[l].length; o++ ) {
						newDependencyParameters[l][o] = new double[context][a];
						if( o != 0 ) {
							dist=Math.min( l-o+1, distance );
						}
						newAncestorMixtureParameters[l][o] = new double[dist];
						context*=a;
					}
				}
			}
			
			dependencyParameters = newDependencyParameters;
			componentMixtureParameters = newComponentMixtureParameters;
			ancestorMixtureParameters = newAncestorMixtureParameters;
			
			if( modLength ) {
				init();
			} else {
				precompute();
			}
		}
		return true;
	}
	
	public void fillInfixScore(int[] seq, int start, int length, double[] scores) {
		
		for( int l = start; l < start+length; l++ ) {
			int current = seq[l];
			localMixtureScore[0] = componentMixtureParameters[l][0] - componentMixtureLogNorm[l] + dependencyParameters[l][0][0][current] - dependencyLogNorm[l][0][0];
			for( int c = 1; c < componentMixtureParameters[l].length; c++ ) {
				int k=ancestorMixtureParameters[l][c].length;
				int an=getOffset(seq, l, c, dependencyParameters[0][0][0].length);
				for( int m = 0; m < k; m++ ) {
					an = next(an, seq, 0, l, c, m );
					ancestorScore[c][m] = ancestorMixtureParameters[l][c][m] - ancestorMixtureLogNorm[l][c] + dependencyParameters[l][c][an][current] - dependencyLogNorm[l][c][an]; 
				}
				localMixtureScore[c] = componentMixtureParameters[l][c] - componentMixtureLogNorm[l] + Normalisation.getLogSum(0, k, ancestorScore[c]);
			}
			scores[l] = Normalisation.getLogSum(0, componentMixtureParameters[l].length, localMixtureScore );
		}
		
	}
	
}