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
package projects.dispom;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.AbstractMap.SimpleEntry;
import java.util.Map.Entry;

import projects.dispom.PFMComparator.NormalizedEuclideanDistance;
import projects.dispom.PFMComparator.PFMDistance;
import de.jstacs.DataType;
import de.jstacs.WrongAlphabetException;
import de.jstacs.algorithms.optimization.ConstantStartDistance;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.classifier.differentiableSequenceScoreBased.OptimizableFunction;
import de.jstacs.classifier.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifier.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifier.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifier.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifier.differentiableSequenceScoreBased.gendismix.LogGenDisMixFunction;
import de.jstacs.classifier.differentiableSequenceScoreBased.logPrior.CompositeLogPrior;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.DataSetKMerEnumerator;
import de.jstacs.data.DiscreteSequenceEnumerator;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.RecyclableSequenceEnumerator;
import de.jstacs.data.Sequence;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.PermutedSequence;
import de.jstacs.data.sequences.annotation.MotifAnnotation;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.FileManager;
import de.jstacs.io.InfixStringExtractor;
import de.jstacs.io.LimitedStringExtractor;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.io.XMLParser;
import de.jstacs.motifDiscovery.KMereStatistic;
import de.jstacs.motifDiscovery.MotifDiscoverer;
import de.jstacs.motifDiscovery.MutableMotifDiscoverer;
import de.jstacs.motifDiscovery.MutableMotifDiscovererToolbox;
import de.jstacs.motifDiscovery.SignificantMotifOccurrencesFinder;
import de.jstacs.motifDiscovery.MutableMotifDiscovererToolbox.InitMethodForDiffSM;
import de.jstacs.motifDiscovery.SignificantMotifOccurrencesFinder.RandomSeqType;
import de.jstacs.motifDiscovery.history.CappedHistory;
import de.jstacs.motifDiscovery.history.NoRevertHistory;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetTagger;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.NormalizedDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.MarkovModelDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.InhomogeneousMarkov;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.UniformHomogeneousDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mix.StrandDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mix.StrandDiffSM.InitMethod;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mix.motif.DurationDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mix.motif.ExtendedZOOPSDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mix.motif.MixtureDurationDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mix.motif.SkewNormalLikeDurationDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mix.motif.UniformDurationDiffSM;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.SafeOutputStream;

/**
 * Discriminative de-novo position distribution and motif finder.
 * 
 * @author Jens Keilwagen
 */
public class Dispom {
	
	// creates a new homogeneous model	
	private static HomogeneousDiffSM getHomSF( AlphabetContainer con, int order, double ess, int length ) throws Exception {
		if( order >= 0 ) {
			return new HomogeneousMMDiffSM( con, order, ess, length );
		} else {
			return new UniformHomogeneousDiffSM( con, ess );
		}		
	}

	// does a heuristic initialization
	private static void doHeuristic( DataSet data[], double[][] dataWeights, double aprioriMean, int length, boolean bothStrands, int maximalMismatches, int candidates, OptimizableFunction f, double weight, DifferentiableStatisticalModel[] score, SafeOutputStream out, int motifID ) throws Exception {
		out.writeln( "heuristic:" );
		//make a statistic over all k-mers near the a priori mean
		DataSet[] selected;
		int start = (int) Math.max( aprioriMean - 250, 0 );
		int end = (int) Math.min( aprioriMean + 250, data[0].getElementLength() );
		if( start == 0 && end == data[0].getElementLength() ) {
			selected = data;
		} else {
			selected = new DataSet[data.length];
			for( int i = 0; i < data.length; i++ ) {
				selected[i] = data[i].getInfixDataSet( start, end-start );
			}
		}
		Hashtable<Sequence, BitSet[]> hash = KMereStatistic.getKmereSequenceStatistic( length, bothStrands, 1, data );
		out.writeln( "found " + hash.size() + " " + length + "-mers" );
		double nFg = data[0].getNumberOfElements(), nBg = data[1].getNumberOfElements();
		//remove those that occur at most as often in the background as in the foreground (relative frequencies) 
		hash = KMereStatistic.removeBackground( hash, 0, 1, 1d/nFg, 1d/nBg );
		out.writeln( "judged " + hash.size() + " " + length + "-mers as interesting" );
		//merge interesting Strings into "families of k-mers"
		Hashtable<Sequence, BitSet[]> merged = KMereStatistic.merge( hash, maximalMismatches, true );
		ArrayList<ComparableElement<Entry<Sequence,BitSet[]>, Double>> list = new ArrayList<ComparableElement<Entry<Sequence,BitSet[]>,Double>>( merged.size() );
		Iterator<Entry<Sequence,BitSet[]>> it = merged.entrySet().iterator();
		Entry<Sequence,BitSet[]> e;
		BitSet[] b;
		double d1, d2, s;
		while( it.hasNext() ) {
			e = it.next();
			b = e.getValue();
			//score each "family"
			d1 = b[0].cardinality()/nFg - b[1].cardinality()/nBg;
			if( d1 > 0 ) {
				list.add( new ComparableElement<Entry<Sequence,BitSet[]>, Double>( e, d1 ) );
			}
		}
		ComparableElement[] ce = list.toArray( new ComparableElement[0] );
		ComparableElement<Entry<Sequence,BitSet[]>,Double> current;
		Arrays.sort( ce );
		
		Sequence seq, rc, cand;
		ArrayList<Sequence> currentSeqs = new ArrayList<Sequence>(), ignore = new ArrayList<Sequence>();
		DoubleList weights = new DoubleList();
		DataSet[] bestSample = new DataSet[1], currentSample = new DataSet[1];
		double[][] bestWeights = new double[1][], currentWeights = new double[1][];
		double val, bestVal = Double.NEGATIVE_INFINITY;
		int[] classIndex = {0}, motifIndex = {motifID};
		boolean[] adjust = {true};
		MutableMotifDiscoverer[] mmd = { (MutableMotifDiscoverer) score[0] };
		int[] len = { mmd[0].getMotifLength( 0 ) };
		for( int j, anz = 0, idx = ce.length-1; anz < candidates && idx >= 0; idx-- ) {
			current = ce[idx];
			seq = current.getElement().getKey();
			s = 0;
			for( j = 0; j < ignore.size(); j++ ) {
				cand = ignore.get( j );
				rc = cand.reverseComplement();
				d1 = seq.getHammingDistance( cand );
				if( bothStrands ) {
					d2 = seq.getHammingDistance( rc );
					d1 = Math.min( d1, d2 );
				}
				if( d1 <= maximalMismatches+1 ) {
					break;
				}
			}
			if( j == ignore.size() ) {
				//test next k-mer that is not to close to previously checked k-mers
				currentSeqs.clear();
				weights.clear();
				//System.out.println( "==================================================" );
				for( j = 0; j < list.size(); j++ ) {
					e = list.get( j ).getElement();
					cand = e.getKey();
					rc = cand.reverseComplement();
					d1 = seq.getHammingDistance( cand );
					if( bothStrands ) {
						d2 = seq.getHammingDistance( rc );
						if( d1 > d2 ) {
							cand = rc;
							d1 = d2;
						}
					}
					if( d1 <= maximalMismatches ) {
						b = e.getValue();
						weights.add( b[0].cardinality()/nFg - b[1].cardinality()/nBg );
						s += b[0].cardinality()/nFg - b[1].cardinality()/nBg;
						currentSeqs.add( cand );
						//System.out.println( cand + "\t" + ( (b[0].cardinality()/nFg - b[1].cardinality()/nBg) / current.getWeight() ) );
					}
				}
				currentSample[0] = new DataSet( "heuristc sample " + anz, currentSeqs.toArray( new Sequence[0] ) );
				weights.multiply( 0, weights.length(), weight/s );
				currentWeights[0] = weights.toArray();
				MutableMotifDiscovererToolbox.initMotif( 0, classIndex, motifIndex, currentSample, currentWeights, adjust, mmd, len, data, dataWeights );
				f.reset();
				val = f.evaluateFunction( f.getParameters( KindOfParameter.PLUGIN ) );
				out.writeln( anz + "\t" + idx + "\t" + seq + "\t" + current.getWeight() + "\t" + val );
				if( val > bestVal ) {
					bestVal = val;
					bestSample[0] = currentSample[0];
					bestWeights[0] = currentWeights[0];
				}
				ignore.add( seq );
				anz++;
			}
		}
		MutableMotifDiscovererToolbox.initMotif( 0, classIndex, classIndex, bestSample, bestWeights, adjust, mmd, len, data, dataWeights );
		out.writeln( "best=" + bestVal );
	}
	
	private static final String[] PREFIX = {
	        DispomParameterSet.HOME, DispomParameterSet.IGNORE_CHAR, DispomParameterSet.FG,	DispomParameterSet.BG,
	        DispomParameterSet.POSITION_DISTR, DispomParameterSet.MEAN, DispomParameterSet.SD,
	        DispomParameterSet.MOTIFS, DispomParameterSet.LENGTH, DispomParameterSet.FLANKING_ORDER, DispomParameterSet.MOTIF_ORDER,
	        DispomParameterSet.FORWARD_STRAND, DispomParameterSet.INITIALIZE,
	        DispomParameterSet.ADJUST_LENGTH, DispomParameterSet.HEURISTIC,
	        DispomParameterSet.LEARNING_PRINCIPLE_KEY, DispomParameterSet.THREADS,	        
	        DispomParameterSet.STARTS, DispomParameterSet.XML_PATH, DispomParameterSet.P_VALUE
	        }; 
	
	private static DataSet getSample( AlphabetContainer con, String fileName, char ignore ) throws FileNotFoundException, WrongAlphabetException, EmptyDataSetException, WrongLengthException, IOException {
		return new DataSet( con, 
				new LimitedStringExtractor(
					new InfixStringExtractor(//TODO
							new SparseStringExtractor( fileName, ignore ) 
					, 900, 100)
				, 100 )
		);
	}
	
	/**
	 * This is the main of Dispom that starts the program. 
	 * 
	 * @param args the arguments for Dispom. Each argument has the form <code>name=value</code>.
	 * 
	 * @throws Exception if something went wrong.
	 */
	public static void main( String[] args ) throws Exception {
		
		//parse parameters
		ParameterSetTagger params = new ParameterSetTagger( PREFIX, new DispomParameterSet() );
		params.fillParameters( "=", args );
		System.out.println( "parameters:" );
		System.out.println( params );
		System.out.println("_________________________________");
		if( !params.hasDefaultOrIsSet() ) {
			System.out.println( "Some of the required parameters are not specified." );
			System.exit( 1 );
		}
		
		//load data
		AlphabetContainer con = DNAAlphabetContainer.SINGLETON;
		String home = params.getValueFromTag( DispomParameterSet.HOME, String.class );
		char ignore = params.getValueFromTag( DispomParameterSet.IGNORE_CHAR, Character.class );

		LearningPrinciple key = params.getValueFromTag( DispomParameterSet.LEARNING_PRINCIPLE_KEY, LearningPrinciple.class );
		int anz = (key==LearningPrinciple.MAP || key==LearningPrinciple.ML)?1:2;
		if( params.isSet( DispomParameterSet.BG ) ) {
			anz = 2;
		}
		
		DataSet[] data = new DataSet[anz];
		data[0] = getSample( con, home + File.separatorChar + params.getValueFromTag( DispomParameterSet.FG, String.class ), ignore );	
		if( anz > 1 ) {
			if( params.isSet( DispomParameterSet.BG ) ) {
				data[1] = getSample( con, home + File.separatorChar + params.getValueFromTag( DispomParameterSet.BG, String.class ), ignore );				
			} else {
				Sequence[] seqs = data[0].getAllElements();
				for( int n = 0; n < seqs.length; n++ ) {
					seqs[n] = new PermutedSequence( seqs[n] );
				}
				data[1] = new DataSet( "permuted foreground sequences", seqs );
			}
			
			//eliminate intersection
			DataSet s = DataSet.diff( data[1], data[0] );
			if( s.getNumberOfElements() < data[1].getNumberOfElements() ) {
				System.out.println( "removed " + (data[1].getNumberOfElements()-s.getNumberOfElements()) + " sequences from background file" );
				data[1] = s;
			}
		}

		//create simple weights
		double[][] weights = new double[anz][];
		for( int i = 0; i < data.length; i++ ) {
			weights[i] = new double[data[i].getNumberOfElements()];
			Arrays.fill( weights[i], 1 );
			System.out.println( i + "\t# = " + data[i].getNumberOfElements() + "\tlength = " + data[i].getElementLength() + "\t" + data[i].getAnnotation() );
		}
		int sl = data[0].getElementLength();		
		
		// create functions
		int starts = 1,
			flOrder = params.getValueFromTag( DispomParameterSet.FLANKING_ORDER, Integer.class ),
			fgOrder = params.getValueFromTag( DispomParameterSet.MOTIF_ORDER, Integer.class ),
			motifL = params.getValueFromTag( DispomParameterSet.LENGTH, Integer.class ),
			motifs = params.getValueFromTag( DispomParameterSet.MOTIFS, Integer.class );
		double essMotif = 4, essNonMotif = 1, f;
		boolean free = false;

		DifferentiableStatisticalModel[] score = new DifferentiableStatisticalModel[anz];
		int n0 = data[0].getNumberOfElements(), n1 = data[1].getNumberOfElements(); 
		if( n0 >= n1 ) {
			f = Math.round(n0/(double)n1);
		} else {
			f = Math.round(n1/(double)n0);
		}
		if( anz > 1 ) {
			score[1] = getHomSF( con, flOrder, f*(essMotif + motifs*essNonMotif), sl );
		}
		HomogeneousDiffSM flanking = getHomSF( con, flOrder, essNonMotif, sl );
		
		int threads = params.getValueFromTag( DispomParameterSet.THREADS, Integer.class );
		double epsilon = 1E-7, lineps = 1E-10, startD = 1;
		AbstractTerminationCondition eps = new SmallDifferenceOfFunctionEvaluationsCondition( epsilon );
		byte algo = Optimizer.QUASI_NEWTON_BFGS;
		SafeOutputStream stream = SafeOutputStream.getSafeOutputStream( System.out );//TODO
		LogGenDisMixFunction objective, initObjective;
		double[] beta = LearningPrinciple.getBeta( key );
		boolean adjust = params.getValueFromTag( DispomParameterSet.ADJUST_LENGTH, Boolean.class );
		
		// create scoring function with hidden motif
		SkewNormalLikeDurationDiffSM motifPenalty = new SkewNormalLikeDurationDiffSM( 1, 50, -1.5, -2.5, 4.2 );
		//System.out.println( motifPenalty );
		System.out.println();
		DifferentiableStatisticalModel motif = new MarkovModelDiffSM( con, motifL, essMotif, true, new InhomogeneousMarkov(fgOrder), motifPenalty );
		Double forwardStrand = params.getValueFromTag( DispomParameterSet.FORWARD_STRAND, Double.class );
		boolean bothStrands = forwardStrand == null || forwardStrand < 1;
		if( bothStrands ) {
			if( forwardStrand == null ) {
				motif = new StrandDiffSM( motif, 0.5, starts, true, InitMethod.INIT_FORWARD_STRAND );
			} else {
				motif = new StrandDiffSM( motif, starts, true, InitMethod.INIT_FORWARD_STRAND, forwardStrand );
			}
		}
		motif.initializeFunctionRandomly( free );
		motif = new NormalizedDiffSM( motif, 1 );
		
		DurationDiffSM pos = null;
		double seqs = 1, sd = params.getValueFromTag( DispomParameterSet.SD, Double.class ), skewSd = 1, mu = params.getValueFromTag( DispomParameterSet.MEAN, Double.class );
		switch( params.getValueFromTag( DispomParameterSet.POSITION_DISTR, PositionDistribution.class ) ) {
			case UNIFORM:
				pos = new UniformDurationDiffSM( 0, sl-motifL, essMotif );
				break;
			case SKEW_NORMAL:
				pos = new SkewNormalLikeDurationDiffSM( 0, sl-motifL, true, mu, 500, true, seqs/2d, seqs/2d*(sd*sd), true, 0, skewSd, starts);
				break;
			case MIXTURE:
				pos = new MixtureDurationDiffSM( starts,
					new UniformDurationDiffSM( 0, sl-motifL, essMotif-seqs ),
					new SkewNormalLikeDurationDiffSM( 0, sl-motifL, true, mu, 500, true, seqs/2d, seqs/2d*(sd*sd), true, 0, skewSd, starts) );
				break;
		}
		ExtendedZOOPSDiffSM fg;
		if( motifs > 1 ) {
			fg = new ExtendedZOOPSDiffSM( ExtendedZOOPSDiffSM.CONTAINS_SOMETIMES_A_MOTIF, sl, starts, true,	flanking, 
					ArrayHandler.createArrayOf( motif, motifs ), ArrayHandler.createArrayOf( pos, motifs ), true );
					
		} else {		
			fg = new ExtendedZOOPSDiffSM( ExtendedZOOPSDiffSM.CONTAINS_SOMETIMES_A_MOTIF, sl, starts, true,	flanking, motif, pos, true );
		}
		score[0] = fg;
		
		LearningPrinciple initKey = (beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX]>0) ? LearningPrinciple.MCL : LearningPrinciple.ML;	
		String initMethod = params.getValueFromTag( DispomParameterSet.INITIALIZE, String.class );
		int restarts = params.getValueFromTag(DispomParameterSet.STARTS, Integer.class );
		if( restarts > 1 && !initMethod.startsWith( "best-random" ) ) {
			System.out.println( "WARNING: set number of starts to 1, since init-method is \"best-random\"" );
			restarts = 1;
		}
		
		String v = initMethod.substring( initMethod.indexOf( '=' ) + 1 );
		double defaultWeight = 6;
		double[][] res, best = null;
		initObjective = new LogGenDisMixFunction( threads, score, data, weights, new CompositeLogPrior(), LearningPrinciple.getBeta( initKey ), true, free );
		objective = new LogGenDisMixFunction( threads, score, data, weights, new CompositeLogPrior(), beta, true, free );
		
		//repeated starts
		DifferentiableStatisticalModel[] current, bestNSF = null;
		for( int r = 0; r < restarts; r++ ) {
			System.out.println( "start " + r + " -------------------------------------------------" );
			current = ArrayHandler.clone( score );		
			
			objective.reset( current );
			initObjective.reset( current );
			
			// initialize
			if( initMethod.startsWith( "best-random" ) ) {
				InitMethodForDiffSM[] initMeth = {InitMethodForDiffSM.NOTHING, InitMethodForDiffSM.NOTHING};
				if( initMethod.startsWith( "best-random=" ) ) {
					initMeth[0] = InitMethodForDiffSM.RANDOMLY;
				} else if( initMethod.startsWith( "best-random-plugin=" ) ) {
					initMeth[0] = InitMethodForDiffSM.PLUG_IN;
				} else if( initMethod.startsWith( "best-random-motif=" ) ) {
					initMeth[0] = InitMethodForDiffSM.MOTIF_RANDOMLY;
				} 
				
				if( initMeth[0] == InitMethodForDiffSM.NOTHING ) {
					throw new IllegalArgumentException( "Initialization method not correctly set." );
				}
					
				ComparableElement<double[], Double>[] pars = MutableMotifDiscovererToolbox.getSortedInitialParameters( current, initMeth, initObjective, Integer.parseInt( v ), stream, 0 );
				objective.setParams( pars[pars.length-1].getElement() );
			} else if( initMethod.startsWith( "specific=" ) ) {
				DataSet[] spec = new DataSet[1];
				double[][] w = new double[1][];
				if( new File( v ).exists() ) {
					spec[0] = new DataSet( con, new SparseStringExtractor( v, ignore ) );
				} else {
					spec[0] = new DataSet( "one sequence", Sequence.create( con, v.trim() ) );
					w[0] = new double[]{ defaultWeight };
				}
				MutableMotifDiscovererToolbox.initMotif( 0, new int[]{0}, new int[]{0}, spec, w, new boolean[]{true}, new MutableMotifDiscoverer[]{(MutableMotifDiscoverer) current[0]}, new int[]{motifL}, data, weights );
			} else {
				for( int i = 0; i < fg.getNumberOfMotifs(); i++ ) {
					if( initMethod.startsWith( "enum" ) ) {
						RecyclableSequenceEnumerator enumeration;
						int k = Integer.parseInt( v );
						if( initMethod.startsWith( "enum-all=" ) ) {
							enumeration = new DiscreteSequenceEnumerator(con, k, bothStrands );
						} else if( initMethod.startsWith( "enum-data=" ) ) {
							enumeration = new DataSetKMerEnumerator( data[0], k, bothStrands );
						} else{
							throw new IllegalArgumentException( "Initialization method not correctly set." );
						}
						System.out.println( "best seed: " + MutableMotifDiscovererToolbox.enumerate( current, 0, i, enumeration, defaultWeight, initObjective, System.out ) );			
					} else if( initMethod.startsWith( "heuristic" ) ) {				
						doHeuristic( data, weights, mu, (int)(0.7 * motifL), bothStrands, 1, Integer.parseInt( v ), initObjective, defaultWeight, current, stream, i );
					} else{
						throw new IllegalArgumentException( "Initialization method unknown." );
					}
				}
			}
			System.out.println( current[0] );
			System.out.println("_________________________________");
	
			// optimize
			objective.reset( current );
			res = MutableMotifDiscovererToolbox.optimize( current, objective, algo, eps, lineps, new ConstantStartDistance(startD),
					stream, false, 
					new CappedHistory( 15, new NoRevertHistory(true,adjust,adjust) ),
					KindOfParameter.PLUGIN, params.getValueFromTag( DispomParameterSet.HEURISTIC, Boolean.class ) );
			
			// select
			if( best == null || res[0][0] > best[0][0] ) {
				bestNSF = current;
				best = res;
			}
		}
		// save classifier
		GenDisMixClassifierParameterSet cps = new GenDisMixClassifierParameterSet( con, sl, algo, epsilon, lineps, startD, free, KindOfParameter.PLUGIN, true, threads );
		GenDisMixClassifier cl = new GenDisMixClassifier( cps, new CompositeLogPrior(), best[0][0], beta, bestNSF );
		cl.setClassWeights( false, best[1] );
		StringBuffer sb = new StringBuffer( 100000 );
		XMLParser.appendObjectWithTags( sb, cl, "classifier" );
		String fName = params.getValueFromTag( DispomParameterSet.XML_PATH, String.class );
		if( !fName.endsWith( ".xml" ) ) {
			fName = fName + ".xml";
		}
		FileManager.writeFile( new File( fName ), sb );
		initObjective.stopThreads();
		objective.stopThreads();

		// show
		System.out.println("_________________________________");
		System.out.println( bestNSF[0] );
		System.out.println( "result: " + best[0][0] );
		
		// predict BSs
		if( params.isSet( DispomParameterSet.P_VALUE ) ) {
			System.out.println("_________________________________");
			SignificantMotifOccurrencesFinder smof;
			if( anz > 1 ) {
				smof = new SignificantMotifOccurrencesFinder( (MotifDiscoverer) bestNSF[0], data[1], weights[1], params.getValueFromTag( DispomParameterSet.P_VALUE, Double.class ) );			
			} else {
				smof = new SignificantMotifOccurrencesFinder( (MotifDiscoverer) bestNSF[0], RandomSeqType.PERMUTED, true, 1000, params.getValueFromTag( DispomParameterSet.P_VALUE, Double.class ) );
			}
			Sequence seq, site, adjusted;
			MotifAnnotation[] ma;
			System.out.println( "prediction" );
			System.out.println();
			System.out.println( "sequence\tposition\tstrand\tbinding site\tadjusted binding site\tp-value" );
			System.out.println( "------------------------------------------------------------------------");
			LinkedList<Sequence> list = new LinkedList<Sequence>();
			for( int j, i = 0; i < data[0].getNumberOfElements(); i++ ) {
				seq = data[0].getElementAt( i );
				ma = smof.findSignificantMotifOccurrences( 0, seq, 0 );
				if( !(ma == null || ma.length == 0) ) {
					for( j = 0; j < ma.length; j++ ) {
						site = seq.getSubSequence( ma[j].getPosition(), ma[j].getLength() );
						if( ma[j].getStrandedness() == Strand.REVERSE ) {
							adjusted = site.reverseComplement();
						} else {
							adjusted = site;
						}
						System.out.println( i + "\t" + ma[j].getPosition() + "\t" + ma[j].getStrandedness() + "\t" + site + "\t" + adjusted + "\t" + ma[j].getAnnotations()[1].getValue() );
						list.add( adjusted );						
					}
				}
			}
			
			double[][] pfm = null;
			if( list.size() > 0 ) {
				System.out.println();
				System.out.println( "PFM:" );
				pfm = PFMComparator.getPFM( new DataSet( "sites", list.toArray( new Sequence[0] ) ) );
				for( int l = 0; l < pfm.length; l++ ) {
					System.out.print( l );
					for( int a = 0; a < pfm[l].length; a++ ) {
						System.out.print( "\t" + pfm[l][a] );
					}
					System.out.println();
				}
			}
			
			String transfac = "./transfac.dat";
			if( pfm != null && new File( transfac ).exists() ) {
				System.out.println();
				System.out.println( "comparing with TRANSFAC:");
				ArrayList<SimpleEntry<String, double[][]>> library = PFMComparator.readPFMsFromEMBL( transfac, Integer.MAX_VALUE );
				
				PFMDistance dist = new NormalizedEuclideanDistance();
				ComparableElement<String, Double>[] ce = PFMComparator.find( DNAAlphabet.SINGLETON, pfm, library, dist, 7, 2, true, 0.05 );
				for( int i = 0; i < ce.length; i++ ) {
					System.out.println( i + "\t" + ce[i].getWeight() + "\t" +ce[i].getElement() );
				}
			}
		}
	}
}

/**
 * This class is a container for all parameters of Dispom. It also parses the parameter from Strings.
 *  
 * @author Jens Keilwagen
 */
class DispomParameterSet extends ParameterSet {

	public static final String HOME = "home";
	public static final String IGNORE_CHAR = "ignore";
	public static final String FG = "fg";
	public static final String BG = "bg";
	public static final String MEAN = "mean";
	public static final String SD = "sd";
	public static final String LENGTH = "length";
	public static final String FLANKING_ORDER = "flankOrder";
	public static final String MOTIF_ORDER = "motifOrder";
	public static final String FORWARD_STRAND = "forwardStrand";
	public static final String INITIALIZE = "init";
	public static final String XML_PATH = "xml";
	public static final String ADJUST_LENGTH = "adjust";
	public static final String POSITION_DISTR = "position";
	public static final String LEARNING_PRINCIPLE_KEY = "learning";
	public static final String P_VALUE = "p-val";
	public static final String MOTIFS = "motifs";
	public static final String THREADS = "threads";
	public static final String HEURISTIC = "maxPos";
	public static final String STARTS = "starts";
	
	/*
	 * 0 home
	 * 1 ignore
	 * 2 fg file
	 * 3 bg file
	 * 4 a priori TFBS mean
	 * 5 a priori standard deviation
	 * 6 initial length
	 * 7 flanking order
	 * 8 foreground order
	 * 9 both strands
	 * 10 initialization (best-random=,specific=,...)
	 * 11 xml file name
	 * 12 adjust motif length
	 * 13 position distribution
	 * 14 learning principle key
	 * 15 p-val
	 * 16 # motifs
	 * 17 # threads
	 * 18 maxPos
	 * 19 starts
	 */
	public DispomParameterSet() throws Exception {
		super();

		parameters.add( new SimpleParameter( DataType.STRING, "home directory", "the path to the data directory", true, "./" ) );
		parameters.add( new SimpleParameter( DataType.CHAR, "the ignore char for the data files", "the char that is used to mask comment lines in data files, e.g., '>' in a FASTA-file", true, '>' ) );
		parameters.add( new SimpleParameter( DataType.STRING, "foreground file", "the file name of the foreground data file (the file containing sequences which are expected to contain binding sites of a common motif)", true ) );
		parameters.add( new SimpleParameter( DataType.STRING, "background file", "the file name of the background data file", false ) );
		
		parameters.add( new EnumParameter( PositionDistribution.class, "a switch whether to use uniform, skew-normal, or mixture position distribution", true, PositionDistribution.MIXTURE.name() ) );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "mean", "the mean of the a priori TFBS distribution", true, 250d ) );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "sd", "the sd of the a priori TFBS distribution", true, new NumberValidator<Double>(1d,Double.POSITIVE_INFINITY), 150d ) );
		parameters.add( new SimpleParameter( DataType.INT, "number of motifs", "the number of motifs to be searched for", true, new NumberValidator<Integer>(1,5), 1 ) );
		parameters.add( new SimpleParameter( DataType.INT, "initial motif length", "the motif length that is used at the beginning", true, new NumberValidator<Integer>(1,50), 15 ) );
		parameters.add( new SimpleParameter( DataType.INT, "Markov order for flanking models", "the Markov order of the model for the flanking sequence and the background sequence", true, new NumberValidator<Integer>(-1,5), 0 ) );
		parameters.add( new SimpleParameter( DataType.INT, "Markov order for motif model", "the Markov order of the motif model", true, new NumberValidator<Integer>(0,3), 0 ) );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "forward strand", "the probability for finding the binding sites on the forward strand, if not set this probability is inferred from the data", false, new NumberValidator<Double>(0d,1d) ) );
		parameters.add( new SimpleParameter( DataType.STRING, "initialization method", "the method that is used for initialization, one of 'best-random=<number>', 'best-random-plugin=<number>', 'best-random-motif=<number>', 'enum-all=<length>', 'enum-data=<length>', 'heuristic=<number>', and 'specific=<sequence or file of sequences>'", true ) );
		
		parameters.add( new SimpleParameter( DataType.BOOLEAN, "adjust motif length", "a switch whether to adjust the motif length, i.e., either to shrink or expand", true, true ) );
		parameters.add( new SimpleParameter( DataType.BOOLEAN, "max pos heuristic", "a switch whether to use max. pos. in the heuristic or not", true, true ) );
		
		parameters.add( new EnumParameter( LearningPrinciple.class, "a switch for the learning principle", true, LearningPrinciple.MSP.name() ) );
		parameters.add( new SimpleParameter( DataType.INT, "compute threads", "the number of threads that are use to evaluate the objective function and its gradient", true, new NumberValidator<Integer>(1,128), 4 ) );
		parameters.add( new SimpleParameter( DataType.INT, "starts of Dispom", "the number of independent starts of Dispom", true, new NumberValidator<Integer>(1,100), 1 ) );
		parameters.add( new SimpleParameter( DataType.STRING, "classifier xml-file", "the file name of the xml file containing the classifier", true, "./classifier.xml" ) );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "p-value", "a p-value for predicting binding sites", false, new NumberValidator<Double>(0d,1d) ) );
	}
}

/**
 * This enum indicates which distribution will be used for the position.
 * 
 * @author Jens Keilwagen
 */
enum PositionDistribution {	
	/**
	 * Uniform
	 */
	UNIFORM,
	/**
	 * Skew-normal
	 */
	SKEW_NORMAL,
	/**
	 * Mixture
	 */
	MIXTURE;
}


