package de.jstacs.tools;

import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.tools.ui.galaxy.Galaxy;

/**
 * Interface for the protocol of the run of a {@link JstacsTool}.
 * Within a {@link JstacsTool}, only this generic interface should be used.
 * The implementation of the methods highly depends on the interface, e.g., {@link CLI} or {@link Galaxy}.
 * 
 * @author Jan Grau
 *
 */
public interface Protocol {

	/**
	 * Appends a general message to the protocol
	 * @param str the message
	 */
	public void append(String str);
	
	/**
	 * Appends a heading to the protocol, which is generally highlighted in some appropriate way.
	 * @param heading the heading
	 */
	public void appendHeading(String heading);
	
	/**
	 * Appends a warning to the protocol, which is generally highlighted in some appropriate way.
	 * @param warning the warning
	 */
	public void appendWarning(String warning);
	
	/**
	 * Appends a {@link Throwable} to the protocol. The message of the protocol is generally highlighted in some appropriate way.
	 * Depending on the implementation, the stack trace of the {@link Throwable} may also be displayed.
	 * @param th the {@link Throwable}
	 * @see Throwable#printStackTrace()
	 */
	public void appendThrowable(Throwable th);
	
	/**
	 * Appends some verbatim text (i.e., text that is displayed "as is" regardless of default formatting) to the protocol.
	 * @param verbatim the verbatim text
	 */
	public void appendVerbatim(String verbatim);
}
