package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class BismarkPreProcessingPseudoCounts {//calutaion of Pseudocounts on a-priori-probability Of methylated C in Genome

	public static void main(String[] args) throws IOException {
		String bismarkFile=args[0];
		
		BufferedReader BR=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(bismarkFile)))));
		String line="";
		double sumMethylated=0;
		double sumUnmethylated=0;
		double percMethylatedCInGenome=0.0;
		while ((line = BR.readLine()) != null){
			//<chromosome>	<start position>	<end position>	<methylation percentage>	<count methylated>	<count unmethylated>
			//chromosome02    359     359     83.3333333333333        5       1
			//CM/(CM+CU)
			
				String[] splitLine=line.split("\t");
				sumMethylated+=Integer.parseInt(splitLine[4]);
				sumUnmethylated+=Integer.parseInt(splitLine[5]);
				percMethylatedCInGenome=sumMethylated/(sumMethylated+sumUnmethylated);
		}
		BR.close();
		System.out.println("sumMethylated: "+sumMethylated);
		System.out.println("sumUnmethylated: "+sumUnmethylated);
		System.out.println("percMethylatedCInGenome: "+percMethylatedCInGenome);
		
	}

}
