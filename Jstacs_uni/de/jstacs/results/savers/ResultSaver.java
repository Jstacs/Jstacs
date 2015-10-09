package de.jstacs.results.savers;

import java.io.File;

import de.jstacs.results.Result;


public interface ResultSaver<T extends Result> {

	public boolean isAtomic();
	
	public String[] getFileExtensions(T result);
	
	public boolean writeOutput(T result, File path);
	
	public boolean writeOutput(T result, StringBuffer buf);
	
}
