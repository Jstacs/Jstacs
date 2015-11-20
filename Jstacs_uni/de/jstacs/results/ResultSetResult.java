package de.jstacs.results;

import de.jstacs.io.NonParsableException;

/**
 * Class for a {@link Result} that contains a single {@link ResultSet}. Useful for building hierarchies of {@link Result}s and
 * {@link ResultSet}s.
 * 
 * @author Jan Grau
 *
 */
public class ResultSetResult extends ListResult {

	/**
	 * Creates a new {@link ResultSetResult} with given name, comment, annotation, and content.
	 * @param name the name of the result
	 * @param comment a comment on the result
	 * @param annotation annotation (may be <code>null</code>)
	 * @param result the result
	 */
	public ResultSetResult( String name, String comment, ResultSet annotation, ResultSet result ) {
		super( name, comment, annotation, result );
	}

	/**
	 * Creates a {@link ResultSetResult} from its XML representation.
	 * @param representation the XML representation
	 * @throws NonParsableException if XML could not be parsed
	 */
	public ResultSetResult( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

}
