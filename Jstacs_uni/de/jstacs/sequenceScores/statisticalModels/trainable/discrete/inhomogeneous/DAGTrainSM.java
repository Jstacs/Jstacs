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

package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.data.sequences.ByteSequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.XMLParser;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.ConstraintManager;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.IDGTrainSMParameterSet;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * The abstract class for <b>d</b>irected <b>a</b>cyclic <b>g</b>raphical models
 * ({@link DAGTrainSM}).
 * 
 * @author Jens Keilwagen
 */
public abstract class DAGTrainSM extends InhomogeneousDGTrainSM {

	/**
	 * The constraints for the model.
	 */
	protected InhCondProb[] constraints;

	/**
	 * This is the main constructor. It creates a new {@link DAGTrainSM} from the
	 * given {@link IDGTrainSMParameterSet}.
	 * 
	 * @param params
	 *            the given parameter set
	 * 
	 * @throws CloneNotSupportedException
	 *             if the parameter set could not be cloned
	 * @throws IllegalArgumentException
	 *             if the parameter set is not instantiated
	 * @throws NonParsableException
	 *             if the parameter set is not parsable
	 * 
	 * @see InhomogeneousDGTrainSM#InhomogeneousDGTrainSM(IDGTrainSMParameterSet)
	 */
	protected DAGTrainSM( IDGTrainSMParameterSet params ) throws CloneNotSupportedException, IllegalArgumentException, NonParsableException {
		super( params );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link DAGTrainSM} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link DAGTrainSM} could not be reconstructed out of the
	 *             XML representation (the {@link StringBuffer} could not be
	 *             parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see InhomogeneousDGTrainSM#InhomogeneousDGTrainSM(StringBuffer)
	 */
	protected DAGTrainSM( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.InhomogeneousDGTrainSM#clone()
	 */
	@Override
	public DAGTrainSM clone() throws CloneNotSupportedException {
		DAGTrainSM clone = (DAGTrainSM)super.clone();
		if( constraints != null ) {
			clone.constraints = ArrayHandler.clone( constraints );
		}
		return clone;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#emitDataSet(int, int[])
	 */
	@Override
	public DataSet emitDataSet( int n, int... lengths ) throws NotTrainedException, Exception {
		if( !( lengths == null || lengths.length == 0 ) ) {
			throw new Exception( "This is an inhomogeneous model. Please check parameter lengths." );
		}
		if( !trained ) {
			throw new NotTrainedException();
		}
		// topsort on the graph
		boolean[] visited = new boolean[length];
		int[] order = new int[length];
		int counter1 = 0, counter2, index = 0, l;
		Arrays.fill( visited, false );
		for( counter1 = 0; counter1 < length; counter1++ ) {
			if( constraints[counter1].getMarginalOrder() == 1 ) {
				visited[counter1] = true;
				order[index++] = counter1;
			}
		}

		while( index < length ) {
			for( counter1 = 0; counter1 < length; counter1++ ) {
				if( !visited[counter1] ) {
					l = constraints[counter1].getMarginalOrder() - 1;
					counter2 = 0;
					while( counter2 < l && visited[constraints[counter1].getPosition( counter2 )] ) {
						counter2++;
					}
					if( counter2 == l ) {
						order[index++] = counter1;
						visited[counter1] = true;
					}
				}
			}
		}

		// System.out.println( Arrays.toString( order ) );

		byte[] content = new byte[length];
		Sequence[] s = new Sequence[n];
		Random r = new Random();
		for( counter1 = 0; counter1 < n; counter1++ ) {
			//Arrays.fill( content, (byte) -1 );
			for( counter2 = 0; counter2 < length; counter2++ ) {
				// fill content;
				constraints[order[counter2]].getOutput( content, r.nextDouble() );
			}
			s[counter1] = new ByteSequence( alphabets, content );
		}
		return new DataSet( "sampled from " + getInstanceName(), s );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#getLogPriorTerm()
	 */
	public double getLogPriorTerm() throws Exception {
		double pot, p = 0, ess = getESS();
		if( ess > 0 ) {
			if( !trained ) {
				throw new NotTrainedException();
			}
			double anz2;
			int counter1 = 0, counter2, anz1;
			for( ; counter1 < length; counter1++ ) {
				anz1 = constraints[counter1].getNumberOfSpecificConstraints();
				pot = ess / (double)anz1;
				anz2 = alphabetLength[constraints[counter1].getPosition( constraints[counter1].getMarginalOrder() - 1 )];
				p += anz1 / anz2 * Gamma.logOfGamma( anz2 * pot ) - anz1 * Gamma.logOfGamma( pot );
				for( counter2 = 0; counter2 < anz1; counter2++ ) {
					p += pot * constraints[counter1].getLnFreq( counter2 );
				}
			}
			return p;
		} else {
			return 0;
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#getLogProbFor(de.jstacs.data.Sequence, int, int)
	 */
	@Override
	public double getLogProbFor( Sequence sequence, int startpos, int endpos ) throws NotTrainedException, Exception {
		check( sequence, startpos, endpos );
		double erg = constraints[0].getLnFreq( sequence, startpos );
		for( int i = 1; i < length; i++ ) {
			erg += constraints[i].getLnFreq( sequence, startpos );
		}
		return erg;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#getNumericalCharacteristics()
	 */
	public NumericalResultSet getNumericalCharacteristics() {
		return null;
		/*
		if( !trained )
		{
			throw new NotTrainedException();
		}
		double[] sum = { constraints.length, 0 };
		for( int counter1 = 0; counter1 < constraints.length; counter1++ )
		{
			sum[1] += constraints[counter1].getNumberOfSpecificConstraints();
		}
		return new NumericalResultSet( new NumericalResult[]{
				new NumericalResult( "number of constraints", "the number of used cliques", sum[0] ),
				new NumericalResult( "number of specific constraints", "the number of used parameters", sum[1] ) } );
		*/
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.InhomogeneousDGTrainSM#getStructure()
	 */
	@Override
	public String getStructure() throws NotTrainedException {
		if( trained ) {
			StringBuffer all = new StringBuffer( 500 );
			for( int counter1 = 0; counter1 < constraints.length; counter1++ ) {
				all.append( constraints[counter1].toString() + "\n" );
			}
			
			return all.toString();
		} else {
			throw new NotTrainedException();
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.DiscreteGraphicalTrainSM#toString()
	 */
	@Override
	public String toString() {
		String erg = "description: " + getDescription();
		if( trained ) {
			try {
				erg += "\n\nstructure:\n" + getStructure();
				StringBuffer all = new StringBuffer();
				for( int counter1 = 0; counter1 < constraints.length; counter1++ ) {
					all.append( constraints[counter1].getFreqInfo( alphabets ) + "\n" );
				}
				erg += "\nprobabilities:\n" + all.toString();				
			} catch ( NotTrainedException impossible ) {
				System.exit( 1 );
			}
		}
		return erg;
	}

	/**
	 * This method checks whether a given graph is acyclic.
	 * 
	 * @param length
	 *            the sequence length (which corresponds to the number of nodes
	 *            in the graph)
	 * @param graph
	 *            the specified graph
	 * 
	 * @return <code>true</code> if the given graph is acyclic,
	 *         <code>false</code> otherwise
	 */
	protected static boolean checkAcyclic( int length, int[][] graph ) {
		ArrayList<int[]> edges = new ArrayList<int[]>( length );
		int counter1 = 0, counter2;
		boolean[] used = new boolean[length];
		boolean changed = false;
		Arrays.fill( used, false );
		while( counter1 < length ) {
			if( graph[counter1].length > 1 ) {
				edges.add( graph[counter1++] );
			} else {
				used[counter1++] = true;
				changed = true;
			}
		}
		int[] help;
		while( edges.size() > 0 && changed ) {
			changed = false;
			for( counter1 = 0; counter1 < edges.size(); ) {
				help = edges.get( counter1 );
				for( counter2 = 0; counter2 < help.length - 1 && used[help[counter2]]; counter2++ );

				if( counter2 == help.length - 1 ) {
					used[help[help.length - 1]] = true;
					changed = true;
					edges.remove( counter1 );
				} else {
					counter1++;
				}
			}
		}
		return edges.size() == 0;
	}

	/**
	 * This method creates the constraints for a given structure.
	 * 
	 * @param structure
	 *            the specified structure
	 */
	protected void createConstraints( int[][] structure ) {
		constraints = new InhCondProb[length];
		for( int i = 0; i < length; i++ ) {
			constraints[i] = new InhCondProb( structure[i], alphabetLength, structure[i].length > 1 );
		}
	}

	/**
	 * This method draws the parameter of the model from the likelihood or the
	 * posterior, respectively.
	 * 
	 * @param data
	 *            the given data
	 * @param weights
	 *            the weights for the sequences in the data
	 * 
	 * @throws Exception
	 *             if something went wrong while counting or drawing
	 * 
	 * @see ConstraintManager#countInhomogeneous(de.jstacs.data.AlphabetContainer,
	 *      int, DataSet, double[], boolean,
	 *      de.jstacs.sequenceScores.statisticalModels.trainable.discrete.Constraint...)
	 * @see ConstraintManager#drawFreqs(double, InhCondProb...)
	 */
	protected void drawParameters( DataSet data, double[] weights ) throws Exception {
		if( data != null ) {
			ConstraintManager.countInhomogeneous( alphabets, length, data, weights, true, constraints );
		}
		ConstraintManager.drawFreqs( getESS(), constraints );
		trained = true;
	}

	/**
	 * This method estimates the parameter of the model from the likelihood or
	 * the posterior, respectively.
	 * 
	 * @param data
	 *            the data
	 * @param weights
	 *            the weights for the sequences in the data
	 * 
	 * @throws Exception
	 *             if something went wrong while counting or estimating
	 * 
	 * @see DAGTrainSM#drawParameters(DataSet, double[])
	 */
	protected void estimateParameters( DataSet data, double[] weights ) throws Exception {
		ConstraintManager.countInhomogeneous( alphabets, length, data, weights, true, constraints );
		ConstraintManager.computeFreqs( getESS(), constraints );
		trained = true;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.DiscreteGraphicalTrainSM#getFurtherModelInfos()
	 */
	@Override
	protected StringBuffer getFurtherModelInfos() {
		if( trained ) {
			StringBuffer xml = new StringBuffer( 10000 );
			XMLParser.appendObjectWithTags( xml, constraints, "conditionalProb" );
			return xml;
		} else {
			return null;
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.DiscreteGraphicalTrainSM#setFurtherModelInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void setFurtherModelInfos( StringBuffer xml ) throws NonParsableException {
		if( trained ) {
			constraints = XMLParser.extractObjectForTags( xml, "conditionalProb", InhCondProb[].class );
		}
	}
}
