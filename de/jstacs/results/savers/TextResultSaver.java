package de.jstacs.results.savers;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.results.TextResult;


public class TextResultSaver implements ResultSaver<TextResult> {

	public static void register(){
		ResultSaverLibrary.register( TextResult.class, new TextResultSaver() );
	}

	private TextResultSaver() {
	}



	@Override
	public String[] getFileExtensions( TextResult result ) {
		String[] mimes = result.getMime().split( "\\," );
		return mimes;
	}

	@Override
	public boolean writeOutput( TextResult result, File path ) {
		try{
			FileRepresentation rep = result.getValue();
			PrintWriter wr = new PrintWriter( path );
			wr.println( rep.getContent() );
			wr.close();
			return true;
		}catch(IOException e){
			e.printStackTrace();
			return false;
		}
	}

	@Override
	public boolean isAtomic() {
		return true;
	}

	@Override
	public boolean writeOutput( TextResult result, StringBuffer buf ) {
		buf.append( result.getValue().getContent() );
		return true;
	}

}
