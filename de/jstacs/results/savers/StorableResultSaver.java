package de.jstacs.results.savers;

import java.io.File;
import java.io.IOException;

import de.jstacs.io.FileManager;
import de.jstacs.results.StorableResult;

/**
 * Implements a {@link ResultSaver} for {@link StorableResult}. {@link StorableResult} are stored to files as their XML representation
 * using {@link de.jstacs.Storable#toXML()} of the contained {@link de.jstacs.Storable} object. This object is obtained from the {@link StorableResult#getResultInstance()}
 * method.
 * 
 * @author Jan Grau
 */
public class StorableResultSaver implements ResultSaver<StorableResult> {


	/**
	 * Registers this {@link ResultSaver} in the {@link ResultSaverLibrary}
	 */
	public static void register(){
		ResultSaverLibrary.register( StorableResult.class, new StorableResultSaver() );
	}
	
	private StorableResultSaver() { }

	@Override
	public boolean isAtomic() {
		return true;
	}

	@Override
	public String[] getFileExtensions(StorableResult result) {
		return new String[]{"xml"};
	}

	@Override
	public boolean writeOutput(StorableResult result, File path) {
		try{
			FileManager.writeFile(path, result.getResultInstance().toXML());
			return true;
		}catch(IOException e){
			e.printStackTrace();
			return false;
		}
	}

	@Override
	public boolean writeOutput(StorableResult result, StringBuffer buf) {
		buf.append(result.getResultInstance().toXML());
		return true;
	}
	
	

	
	
}
