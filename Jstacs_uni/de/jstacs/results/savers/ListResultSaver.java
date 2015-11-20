package de.jstacs.results.savers;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;

import de.jstacs.results.ListResult;

/**
 * {@link ResultSaver} for {@link ListResult}s.
 * Contents are saved to a text file or {@link StringBuffer} using the {@link ListResult#print(PrintWriter)} method.
 * 
 * @author Jan Grau
 *
 */
public class ListResultSaver implements ResultSaver<ListResult> {

	/**
	 * Registers this {@link ResultSaver} in the {@link ResultSaverLibrary}
	 */
	public static void register(){
		ResultSaverLibrary.register( ListResult.class, new ListResultSaver() );
	}
	
	private ListResultSaver(){}
	
	@Override
	public boolean isAtomic() {
		return true;
	}

	@Override
	public String[] getFileExtensions( ListResult result ) {
		return new String[]{"tsv"};
	}

	@Override
	public boolean writeOutput( ListResult result, File path ) {
		try{
			PrintWriter wr = new PrintWriter( path );
			result.print( wr );
			wr.close();
			return true;
		}catch(IOException e){
			e.printStackTrace();
			return false;
		}
	}

	@Override
	public boolean writeOutput( ListResult result, StringBuffer buf ) {
		StringWriter sw = new StringWriter();
		result.print( new PrintWriter( sw ) );
		buf.append(sw.toString());
		return true;
	}

	
	
	
	
}
