package de.jstacs.results.savers;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import de.jstacs.io.FileManager;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.results.TextResult;

/**
 * {@link ResultSaver} for a {@link TextResult}.
 * The contents of the {@link TextResult} are saved to disk using the {@link FileRepresentation#getContent()} method.
 * If the file defined by the contained {@link FileRepresentation} already exists, it is simply copied to the new location.
 * 
 * 
 * @author Jan Grau
 *
 */
public class TextResultSaver implements ResultSaver<TextResult> {

	/**
	 * Registers this {@link ResultSaver} in the {@link ResultSaverLibrary}
	 */
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
			if(rep.getFilename() != null && (new File(rep.getFilename())).exists()){
				FileManager.copy(rep.getFilename(), path.getAbsolutePath());
			}else{
				PrintWriter wr = new PrintWriter( path );
				wr.println( rep.getContent() );
				wr.close();
			}
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
