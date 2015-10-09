package de.jstacs.results;

import de.jstacs.io.NonParsableException;


public class ResultSetResult extends ListResult {

	public ResultSetResult( String name, String comment, ResultSet annotation, ResultSet result ) {
		super( name, comment, annotation, result );
	}

	public ResultSetResult( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

}
