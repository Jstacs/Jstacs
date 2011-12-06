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

package de.jstacs.models.discrete.inhomogeneous;

import java.util.Map;
import java.util.StringTokenizer;
import java.util.TreeMap;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.data.DataSet;
import de.jstacs.io.XMLParser;
import de.jstacs.models.Model;
import de.jstacs.models.discrete.DGMParameterSet;
import de.jstacs.models.discrete.inhomogeneous.parameters.FSDAGMParameterSet;

/**
 * This class can be used for any discrete <b>f</b>ixed <b>s</b>tructure
 * <b>d</b>irected <b>a</b>cyclic <b>g</b>raphical model ( {@link FSDAGModel}).
 * 
 * @author Jens Keilwagen
 */
public class FSDAGModel extends DAGModel {

	/**
	 * This is the main constructor. It creates a new {@link FSDAGModel} from
	 * the given {@link FSDAGMParameterSet}.
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
	 * @see DAGModel#DAGModel(de.jstacs.models.discrete.inhomogeneous.parameters.IDGMParameterSet)
	 */
	public FSDAGModel( FSDAGMParameterSet params ) throws CloneNotSupportedException, IllegalArgumentException, NonParsableException {
		super( params );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link FSDAGModel} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link FSDAGModel} could not be reconstructed out of
	 *             the XML representation (the {@link StringBuffer} could not be
	 *             parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see DAGModel#DAGModel(StringBuffer)
	 */
	public FSDAGModel( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.Model#getInstanceName()
	 */
	public String getInstanceName() {
		return "fixed structure directed acyclic graphical model";
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.AbstractModel#getMaximalMarkovOrder()
	 */
	@Override
	public byte getMaximalMarkovOrder() {
		int max = constraints[0].getMarginalOrder();
		for( int i = 1; i < length; i++ ) {
			if( max < constraints[i].getMarginalOrder() ) {
				max = constraints[i].getMarginalOrder();
			}
		}
		return (byte)( max - 1 );
	}

	private static final String XML_TAG = "FSDAGModel";

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.DiscreteGraphicalModel#getXMLTag()
	 */
	@Override
	protected String getXMLTag() {
		return XML_TAG;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.Model#train(de.jstacs.data.Sample, double[])
	 */
	public void train( DataSet data, double[] weights ) throws Exception {
		estimateParameters( data, weights );
	}

	/**
	 * Computes the model with structure <code>graph</code>.
	 * 
	 * @param data
	 *            the {@link DataSet}
	 * @param weights
	 *            the weights for the sequences in the {@link DataSet}
	 * @param graph
	 *            the graph
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public void train( DataSet data, double[] weights, int[][] graph ) throws Exception {
		// check
		if( !checkAcyclic( length, graph ) ) {
			throw new IllegalArgumentException( "the graph is not acyclic" );
		}
		trainUnchecked( data, weights, graph );
	}

	/**
	 * This method draws the parameters of the model from the a posteriori
	 * density. For drawing from the prior you have to set the data and their
	 * weights to <code>null</code>. Furthermore this method enables you to
	 * specify a new graph structure.
	 * 
	 * @param data
	 *            a {@link DataSet} or <code>null</code>
	 * @param weights
	 *            the (positive) weights for each sequence of the {@link DataSet}
	 *            or <code>null</code>
	 * @param graph
	 *            the graph or <code>null</code> for the current graph
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see DAGModel#drawParameters(DataSet, double[])
	 * @see DAGModel#checkAcyclic(int, int[][])
	 */
	public void drawParameters( DataSet data, double[] weights, int[][] graph ) throws Exception {
		if( graph != null ) {
			// check
			if( !checkAcyclic( length, graph ) ) {
				throw new IllegalArgumentException( "the graph is not acyclic" );
			}
			params.getParameterAt( 2 ).setValue( FSDAGMParameterSet.encode( graph ) );
			// create
			createConstraints( graph );
		}
		// draw
		drawParameters( data, weights );
	}

	/**
	 * This method trains the model with given graph structure on the given
	 * data. Be careful if you use this method since it does not check if the
	 * graph is acyclic.
	 * 
	 * @param data
	 *            the given data
	 * @param weights
	 *            the (positive) weights for each sequence of data
	 * @param graph
	 *            the given graph
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see DAGModel#checkAcyclic(int, int[][])
	 */
	private void trainUnchecked( DataSet data, double[] weights, int[][] graph ) throws Exception {
		// set parameter
		params.getParameterAt( 2 ).setValue( FSDAGMParameterSet.encode( graph ) );
		// create
		createConstraints( graph );
		// estimate
		estimateParameters( data, weights );
	}

	/**
	 * Computes the models with structure <code>graph</code>.
	 * 
	 * @param models
	 *            an array of {@link de.jstacs.models.AbstractModel}s containing
	 *            only instances of {@link FSDAGModel}
	 * @param data
	 *            the {@link DataSet}
	 * @param weights
	 *            the weights for the sequences in the {@link DataSet}
	 * @param graph
	 *            the graph
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static void train( Model[] models, int[][] graph, double[][] weights, DataSet... data ) throws Exception {
		// check
		if( !checkAcyclic( graph.length, graph ) ) {
			throw new IllegalArgumentException( "the graph is not acyclic" );
		}
		if( data.length == 1 ) {
			for( int i = 0; i < models.length; i++ ) {
				( (FSDAGModel)models[i] ).trainUnchecked( data[0], weights[i], graph );
			}
		} else {
			for( int i = 0; i < models.length; i++ ) {
				( (FSDAGModel)models[i] ).trainUnchecked( data[i], weights[i], graph );
			}
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.inhomogeneous.InhomogeneousDGM#set(de.jstacs.models.discrete.DGMParameterSet, boolean)
	 */
	@Override
	protected void set( DGMParameterSet params, boolean trained ) throws CloneNotSupportedException, NonParsableException {
		super.set( params, trained );
		byte i = 0, j;
		// extract
		StringBuffer help, xml = new StringBuffer( (String)params.getParameterAt( 2 ).getValue() );
		StringTokenizer t;
		int[][] graph = new int[length][];
		Map<String,String> filter = new TreeMap<String, String>();
		for( i = 0; i < length; i++ ) {
			filter.clear();
			filter.put( "node", ""+i );
			help = XMLParser.extractForTag( xml, "parents", null, filter );
			if( help == null ) {
				graph[i] = new int[]{ i };
			} else {
				t = new StringTokenizer( help.toString(), "," );
				graph[i] = new int[t.countTokens() + 1];
				for( j = 0; t.hasMoreTokens(); j++ ) {
					graph[i][j] = Integer.parseInt( t.nextToken().trim() );
					if( graph[i][j] < 0 || graph[i][j] >= length || graph[i][j] == i ) {
						throw new IllegalArgumentException( "The graph was not correct encoded. See parents from node " + i );
					}
				}
				graph[i][j] = i;
				if( XMLParser.extractForTag( xml, "parents", null, filter ) != null ) {
					throw new IllegalArgumentException( "The graph was not correct encoded. There are at least 2 sets of parents for node " + i );
				}
			}
		}
		// check
		if( !checkAcyclic( length, graph ) ) {
			throw new IllegalArgumentException( "the graph is not acyclic" );
		}
		// create
		createConstraints( graph );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.inhomogeneous.DAGModel#getStructure()
	 */
	@Override
	public String getStructure() {
		if( trained ) {
			try {
				return super.getStructure();
			} catch ( NotTrainedException e ) {
				// will never happen
				RuntimeException r = new RuntimeException( e.getMessage() );
				throw r;
			}
		}

		StringBuffer all = new StringBuffer( 500 );
		for( int counter1 = 0; counter1 < constraints.length; counter1++ ) {
			all.append( constraints[counter1].toString() + "\n" );
		}
		return all.toString();
	}
}
