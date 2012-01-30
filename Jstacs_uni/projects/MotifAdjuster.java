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

package projects;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.StringExtractor;
import de.jstacs.motifDiscovery.MotifDiscoverer.KindOfProfile;
import de.jstacs.sequenceScores.statisticalModels.trainable.AbstractTrainableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.homogeneous.HomogeneousMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.homogeneous.parameters.HomMMParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.BayesianNetworkTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.ModelType;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.BayesianNetworkTrainSMParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.StrandTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Parameterization;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.motif.ZOOPSTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.motif.positionprior.GaussianLikePositionPrior;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.motif.positionprior.PositionPrior;

/**
 * This class provides a main that is used for the MotifAdjuster .
 * 
 * @author Jens Keilwagen
 */
public class MotifAdjuster
{
	/**
	 * @param args
	 *            <ol>
	 *            <li> file: the location of the data set (String)
	 *            <li> ignoreChar: char for comment lines (e.g. for a FastA-file '&gt;') (char)
	 *            <li> length: the motif length (int)
	 *            <li> fgOrder: the order of the inhomogeneuous Markov model that is uses for the motif; 0 yields in a
	 *            PWM (byte)
	 *            <li> ess: the equivalent sample size that is used for the mixture model (double &gt;= 0)
	 *            <li> bothStrands: use both strands (boolean)
	 *            <li> output: output of the EM (boolean)
	 *            <li> sigma: the sigma of the truncated discrete Gaussian distribution (double&gt;0)
	 *            <li> p(no motif): the probability for finding no motif (0&lt;=double&lt;1)
	 *            </ol>
	 */
	public static void main( String[] args )
	{
		System.out.println( "java ... MotifAdjuster <file> <ignoreChar> <length> <fgOrder> <ess> <bothStrands> <output> <sigma> <p(no motif)>" );
		try
		{
			AlphabetContainer con =DNAAlphabetContainer.SINGLETON;

			char ignore = args[1].charAt( 0 );
			DataSet s = new DataSet( con, new StringExtractor( new File( args[0] ), 200, ignore ) );

			if( s.getElementLength() == 0 )
			{
				System.out.println( "All sequences have to have the same length." );
			}
			else
			{
				System.out.println( s.getAnnotation() + ": " + s.getNumberOfElements() + " sequences of length "
						+ s.getElementLength() );

				int sl = Integer.parseInt( args[2] );
				int l = s.getElementLength(), max = (l - sl) / 2;

				byte fgOrder = Byte.parseByte( args[3] );
				double ess = Double.parseDouble( args[4] );
				boolean bothStrands = Boolean.parseBoolean( args[5] );
				double sigma = Double.parseDouble( args[7] );
				double pBg = Double.parseDouble( args[8] ), pwmESS = (1d-pBg)*ess, strandHyper=pwmESS/2d;
				
				BayesianNetworkTrainSMParameterSet p = new BayesianNetworkTrainSMParameterSet( con, sl, pwmESS, "foreground model", ModelType.IMM, fgOrder, LearningType.ML_OR_MAP );
				AbstractTrainableStatisticalModel motifModel;
				TerminationCondition stop = new SmallDifferenceOfFunctionEvaluationsCondition( 1E-6 );
				if( bothStrands )
				{
					motifModel = new StrandTrainSM( new BayesianNetworkTrainSM( p ), 1, new double[]{ strandHyper, strandHyper }, 1, stop, Parameterization.LAMBDA );
					((StrandTrainSM) motifModel).setOutputStream( null );
				}
				else
				{
					motifModel = new BayesianNetworkTrainSM( p );
				}

				AbstractTrainableStatisticalModel backgroundModel = new HomogeneousMM( new HomMMParameterSet( con, (s.getElementLength() - (1d-pBg)*sl) * ess, null, (byte) 0 ) );

				// here you can alter the prior
				PositionPrior pr = new GaussianLikePositionPrior( l, max, sigma );
				System.out.println( "prior prop. to: exp( -(" + max + " - l)^2/(2 * " + sigma + "^2) )" );

				
				ZOOPSTrainSM em;
				
				em = new ZOOPSTrainSM( motifModel, backgroundModel, false, 10, 1d-pBg, pr, 1d, stop, Parameterization.LAMBDA );
				//em = new ZOOPSTrainSM( motifModel, backgroundModel, false, 10, new double[]{4,1}, pr, 1d, stop, Parameterization.LAMBDA );

				if( Boolean.parseBoolean( args[6] ) )
				{
					em.setOutputStream( System.out );
				}
				else
				{
					em.setOutputStream( null );
				}
				em.train( s );
				
				System.out.println();
				System.out.println( "models: " );

				System.out.println( em );
				motifModel = (AbstractTrainableStatisticalModel) em.getModel( 0 );

				Sequence seq;
				int start;
				StrandTrainSM strand = null;
				if( bothStrands )
				{
					strand = (StrandTrainSM) motifModel;
				}
				String annot, line;
				System.out.println( "results for " + s.getNumberOfElements() + " sites" );
				System.out.println();
				System.out.println( "\"annotation\"\tsequence\tcontains BS\tpredicted shift\tpredicted strand\tadjusted BS" );

				BufferedReader reader = new BufferedReader( new FileReader( args[0] ) );
				while( (line = reader.readLine()) != null )
				{
					annot = " ";
					while( line.charAt( 0 ) == ignore )
					{
						annot = line;
						line = reader.readLine();
					}
					seq = Sequence.create( con, line );
					System.out.print( "\"" + annot.substring( 1 ).trim() + "\"\t" + seq + "\t" );
					if( em.getIndexOfMaximalComponentFor( seq ) == 1 )
					{
						System.out.print( "0" );
					}
					else
					{
						double[] prof = em.getProfileOfScoresFor( 0, 0, seq, 0, KindOfProfile.UNNORMALIZED_CONDITIONAL );
						start = getIndexOfMax( prof );
						System.out.print( "1\t" + (start - max) + "\t" );
						seq = seq.getSubSequence( start, sl );

						if( bothStrands )
						{
							if( strand.getIndexOfMaximalComponentFor( seq ) == 0 )
							{
								System.out.print( "forward " );
							}
							else
							{
								System.out.print( "rev. compl." );
								seq = seq.reverseComplement();
							}
						}
						else
						{
							System.out.print( "forward " );
						}
						System.out.print( "\t" + seq );
					}
					System.out.println();
				}
			}
		}
		catch( Exception e )
		{
			e.printStackTrace();
		}
	}
	
	private static int getIndexOfMax( double[] array ) {
		int idx = 0, i = 1;
		for( ; i < array.length; i++ ) {
			if( array[idx] < array[i] ) {
				idx = i;
			}
		}
		return idx;
	}
}