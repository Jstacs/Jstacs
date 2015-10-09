package de.jstacs.tools;

import de.jstacs.parameters.ParameterSet;


public interface JstacsTool {

	public ParameterSet getToolParameters();
	
	public ToolResult run(ParameterSet parameters, Protocol protocol, ProgressUpdater progress) throws Exception;

	public String getToolName();
	
	public String getShortName();
	
	public String getDescription();
	
	public String getHelpText();
	
}
