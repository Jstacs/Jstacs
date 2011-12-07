package de.jstacs;

import de.jstacs.parameters.Parameter;
import de.jstacs.results.Result;

/**
 * Superclass for all Jstacs entities that have a name, a comment, and a data type as annotations.
 * Specifically, such entities are {@link Result}s for storing results of some computation together with an annotation on their
 * meaning, and {@link Parameter}s, which allow for annotating external parameters, especially parameters of constructors.
 * 
 * @author Jan Grau
 *
 */
public abstract class AnnotatedEntity implements Storable{

	/**
	 * The name of the result.
	 */
	protected String name;

	/**
	 * The comment for the result.
	 */
	protected String comment;

	/**
	 * The data type of the result.
	 */
	protected DataType datatype;

	
	/**
	 * Returns the data type of the {@link AnnotatedEntity}.
	 * 
	 * @return the data type of the {@link AnnotatedEntity}
	 */
	public DataType getDatatype() {
		return datatype;
	}

	/**
	 * Returns the value of the {@link AnnotatedEntity}.
	 * 
	 * @return the value of the {@link AnnotatedEntity}
	 */
	public abstract Object getValue();

	/**
	 * Returns the name of the {@link AnnotatedEntity}.
	 * 
	 * @return the name of the {@link AnnotatedEntity}
	 */
	public String getName() {
		return name;
	}

	/**
	 * Returns the comment on the {@link AnnotatedEntity}.
	 * 
	 * @return the comment on the {@link AnnotatedEntity}
	 */
	public String getComment() {
		return comment;
	}
	
}
