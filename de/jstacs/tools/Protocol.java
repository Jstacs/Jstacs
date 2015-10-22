package de.jstacs.tools;


public interface Protocol {

	public void append(String str);
	
	public void appendHeading(String heading);
	
	public void appendWarning(String warning);
	
	public void appendThrowable(Throwable th);
	
	public void appendVerbatim(String verbatim);
}
