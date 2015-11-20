package de.jstacs.tools;

import de.jstacs.cli.CLI;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.utils.galaxy.Galaxy;

/**
 * Interface for a generic Jstacs tool.
 * 
 * Implementations of this interface may be used in the {@link CLI}, {@link Galaxy}, and JavaFX-based GUI classes
 * for creating generic interfaces to tools created using Jstacs
 * 
 * @author Jan Grau
 *
 */
public interface JstacsTool {

	/**
	 * Returns the input parameters of this tool. 
	 * The parameters should be empty but may have default values, which are used in all interface variants (command line, Galaxy, JavaFX GUI).
	 * @return the input parameters
	 */
	public ParameterSet getToolParameters();
	
	/**
	 * Runs the tool using the provided (now filled) parameters, which are in structure identical to those returned by {@link #getToolParameters()}. These parameters should only be used for this run and should not affect subsequent runs of the same tool.
	 * {@link Protocol} and {@link ProgressUpdater} may be used for indicating the tool's progress. Depending on the implementation and interface variant, these may be rendered differently (or not at all). This method returns all results of this tool encapsulated in a {@link ToolResult}.
	 * @param parameters the input parameters
	 * @param protocol the protocol
	 * @param progress the progress updater
	 * @return the results of this tool
	 * @throws Exception
	 */
	public ToolResult run(ParameterSet parameters, Protocol protocol, ProgressUpdater progress) throws Exception;

	/**
	 * Returns a descriptive, human readable name for this tool.
	 * @return the name
	 */
	public String getToolName();
	
	/**
	 * Returns a name (preferably short and without spaces) for referring to this tool on the command line.
	 * @return the short name
	 */
	public String getShortName();
	
	/**
	 * Returns a short description (half a sentence) on what this tool does.
	 * @return the description
	 */
	public String getDescription();
	
	/**
	 * Returns a detailed help text for this tool, describing the purpose of the tool, all parameters and results.
	 * May use <href="http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html">reStructureText</a> markup.
	 * @return the help text
	 */
	public String getHelpText();
	
}
