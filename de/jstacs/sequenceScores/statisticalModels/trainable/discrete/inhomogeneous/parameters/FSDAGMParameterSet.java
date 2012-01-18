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

package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters;

import de.jstacs.DataType;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.FSDAGTrainSM;

/**
 * The class for the parameters of a {@link FSDAGTrainSM} (<b>f</b>ixed
 * <b>s</b>tructure <b>d</b>irected <b>a</b>cyclic <b>g</b>raphical
 * <b>m</b>odel).
 * 
 * @author Jens Keilwagen
 */
public class FSDAGMParameterSet extends IDGTrainSMParameterSet {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link FSDAGMParameterSet} out of its XML representation.
	 * 
	 * @param s
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link FSDAGMParameterSet} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer} could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see IDGTrainSMParameterSet#IDGTrainSMParameterSet(StringBuffer)
	 */
	public FSDAGMParameterSet( StringBuffer s ) throws NonParsableException {
		super( s );
	}

	/**
	 * This constructor creates an empty {@link FSDAGMParameterSet} set for a
	 * {@link FSDAGTrainSM}.
	 * 
	 * @see FSDAGMParameterSet#FSDAGMParameterSet(Class)
	 */
	public FSDAGMParameterSet() {
		this( FSDAGTrainSM.class );
	}

	/**
	 * This constructor creates an {@link FSDAGMParameterSet} instance. It sets
	 * the {@link AlphabetContainer}, the length, the ess (<b>e</b>quivalent
	 * <b>s</b>ample <b>s</b>ize) and the model description as well as a
	 * {@link String} describing the graph structure.
	 * 
	 * @param alphabet
	 *            the {@link AlphabetContainer} that is used in the model
	 * @param length
	 *            the length of the model (has to be positive)
	 * @param ess
	 *            the ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize) of the
	 *            model (has to be positive)
	 * @param description
	 *            a short description of the model (for a better handling of the
	 *            object by the user)
	 * @param graph
	 *            the graph description {@link String}, encodes in XML-like
	 *            manner the parents of each node &quot;&lt;parents
	 *            node=&quot;i&quot;&gt;j,k,l&lt;/parents&gt;&quot;
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see FSDAGMParameterSet#encode(int[][])
	 * @see FSDAGMParameterSet#FSDAGMParameterSet(Class, AlphabetContainer, int,
	 *      double, String, String)
	 */
	public FSDAGMParameterSet( AlphabetContainer alphabet, int length, double ess, String description, String graph ) throws Exception {
		this( FSDAGTrainSM.class, alphabet, length, ess, description, graph );
	}

	/**
	 * This the constructor creates an empty {@link FSDAGMParameterSet} from the
	 * class that can be instantiated using this {@link FSDAGMParameterSet}.
	 * 
	 * @param clazz
	 *            the class of the object that will be created with this
	 *            parameter set
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.FSDAGTrainSM
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.FSDAGModelForGibbsSampling
	 * @see IDGTrainSMParameterSet#IDGTrainSMParameterSet(Class)
	 */
	protected FSDAGMParameterSet( Class<? extends FSDAGTrainSM> clazz ) {
		super( clazz );
		try {
			parameters.add( new SimpleParameter( DataType.STRING,
					"graph structure",
					"the graph structure for the model (has to be acyclic)&lt;br&gt;" + "The graph has to be given in a semi-XML-format, i.e. if 1, 2 and 5 are the parents from node 0 enter&lt;br&gt;"
							+ "&lt;parents node=0&gt;1,2,5&lt;/parents&gt;&lt;br&gt;"
							+ "the root nodes do not have to be insert explicitly",
					false ) );
		} catch ( DatatypeNotValidException doesnothappen ) { }
	}

	/**
	 * This constructor creates an {@link FSDAGMParameterSet} instance for the
	 * specified class. It sets the {@link AlphabetContainer}, the length, the
	 * ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize) and the model
	 * description as well as a {@link String} describing the graph structure.
	 * 
	 * @param clazz
	 *            the class of the object that will be created with this
	 *            parameter set
	 * @param alphabet
	 *            the {@link AlphabetContainer} that is used in the model
	 * @param length
	 *            the length of the model (has to be positive)
	 * @param ess
	 *            the ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize) of the
	 *            model (has to be positive)
	 * @param description
	 *            a short description of the model (for a better handling of the
	 *            object by the user)
	 * @param graph
	 *            the graph description {@link String}, encodes in XML-like
	 *            manner the parents of each node &quot;&lt;parents
	 *            node=&quot;i&quot;&gt;j,k,l&lt;/parents&gt;&quot;
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see FSDAGMParameterSet#encode(int[][])
	 * @see IDGTrainSMParameterSet#IDGTrainSMParameterSet(Class, AlphabetContainer, int,
	 *      double, String)
	 */
	protected FSDAGMParameterSet( Class<? extends FSDAGTrainSM> clazz, AlphabetContainer alphabet, int length, double ess,
									String description, String graph ) throws Exception {
		super( clazz, alphabet, length, ess, description );
		parameters.add( new SimpleParameter( DataType.STRING,
				"graph structure",
				"the graph structure for the model (has to be acyclic)&lt;br&gt;" + "The graph has to be given in a semi-XML-format, i.e. if 1, 2 and 5 are the parents from node 0 enter&lt;br&gt;"
						+ "&lt;parents node=0&gt;1,2,5&lt;/parents&gt;&lt;br&gt;"
						+ "the root nodes do not have to be insert explicitly",
				false ) );
		parameters.get( 2 ).setValue( graph );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
	 */
	@Override
	public String getInstanceComment() {
		return "a directed acyclic graphical (DAG) model with fixed structure (FS)";
	}

	/**
	 * This method can be used to encode an adjacency list to a graph
	 * description {@link String} (e.g. for the different constructors which
	 * requires graph description {@link String}s).
	 * 
	 * @param graph
	 *            <code>graph[i]</code> contains the parents of node
	 *            <code>i</code>
	 * 
	 * @return the graph description {@link String}
	 * 
	 * @see FSDAGMParameterSet#FSDAGMParameterSet(AlphabetContainer, int,
	 *      double, String, String)
	 * @see FSDAGMParameterSet#FSDAGMParameterSet(Class, AlphabetContainer, int,
	 *      double, String, String)
	 */
	public static String encode( int[][] graph ) {
		StringBuffer help, encoded = new StringBuffer( 1000 );
		int i = 0, j;
		for( ; i < graph.length; i++ ) {
			if( graph[i] != null && graph[i].length > 0 ) {
				help = new StringBuffer( 100 );
				help.append( graph[i][0] );
				for( j = 1; j < graph[i].length; j++ ) {
					help.append( "," + graph[i][j] );
				}
				XMLParser.addTagsAndAttributes( help, "parents", "node=\"" + i + "\"" );
				encoded.append( help );
			}
		}
		return encoded.toString();
	}

}
