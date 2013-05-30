package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameter;

/**
 * This class enables you to input your own structure defined by some constraints.
 * 
 * @author Jens Keilwagen
 */
public final class ConstraintParameterSet extends ParameterSet
{
	/**
	 * The main constructor. 
	 */
	public ConstraintParameterSet()
	{
		super();
		initParameterList( 1 );
		// TODO better comments, Validator?
		try {
			parameters
					.add( new SimpleParameter(
							DataType.STRING,
							"the used constraints",
							"There are two ways entering the constraints."
									+ "<ol>"
									+ "<li> Enter <b>single constraints</b> in the following way: node_1, node_2, ..., node_k.</li>"
									+ "<li> Enter <b>groups of constraints</b> in the following way: &quot;m&lt;number of nodes&gt;s&lt;number of skiped nodes&gt;&quot;."
									+ "E.g. m2s1 is the set {0,2; 1,3; 2,4; ...}."//<br>"
									+ "The &quot;x&quot; can be used to determine that each skip should be used (e.g. m2sx).</li>"
									+ "</ol>" + "To separte the constraints use &quot;;&quot;."//<br><br>"
									+ "Additionally for each node a constraint is assumed (i.e. {0; 1; 2; ...}).", true ) );
		} catch ( Exception e ) {
			throw new RuntimeException( e.getCause() );
		}
	}

	/**
	 * The constructor for the <code>Storable</code> interface.
	 * 
	 * @param xml
	 *            the StringBuffer
	 * 
	 * @throws NonParsableException
	 *             if the StringBuffer is not parsable
	 */
	public ConstraintParameterSet( StringBuffer xml ) throws NonParsableException
	{
		super( xml );
	}
}