package projects.talecorrect;

import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ui.cli.CLI;

public class TALEcorrection {

	public static void main( String[] args ) throws Exception {
		
		JstacsTool[] tools = new JstacsTool[]{
				new CorrectTALESequences(),
				new PolishTALESubstitutions()};

		CLI cli = new CLI( tools );
		
		cli.run( args );

	}

}
