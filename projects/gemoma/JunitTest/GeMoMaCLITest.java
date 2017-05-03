package projects.gemoma.JunitTest;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.ArrayList;
import java.util.Date;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import projects.gemoma.ExtractRNAseqEvidence.Stranded;
import projects.gemoma.GeMoMa;
import projects.gemoma.Tools.Ambiguity;

/**
 * 
 * @author Jens Keilwagen
 */
public class GeMoMaCLITest {

	static int start, end;
	
//for creating infrastructure, clean up, and nice things
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		start = end = 0;
		File dir = new File( out + "-"+ System.currentTimeMillis() + File.separator );
		dir.mkdirs();
		out = dir.getAbsolutePath() + File.separator;
		System.out.println(out);
		System.out.println();
		System.out.println("start: " + new Date());
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
		if( end == start ) {
			//if successful then delete
			System.out.println("all tests successful -> delete all files");
			File dir = new File( out );
			for(File file: dir.listFiles()) { 
			        file.delete();
			}
			dir.deleteOnExit();
		}
		System.out.println("end: " + new Date());
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
		System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	}

//own test implementation	
	static String in = "projects/gemoma/JunitTest/given/";
	static String out = "projects/gemoma/JunitTest/results";

	static void cliCheck( ArrayList<String> args, String[] given, String[] newResults ) throws Exception {
		start++;
		assertTrue( "Different number of results", given.length == newResults.length );
		args.add(0, "CLI");
		args.add("outdir="+out);
		GeMoMa.main( args.toArray(new String[args.size()]) );
		
		for( int i = 0; i < given.length; i++ ) {
			GeMoMaTest.assertFile( null, in + given[i], out+newResults[i] );
		}
		end++;
	}
	
	@Test
	public void testExtractor() throws Exception {
		ArrayList<String> list = new ArrayList<String>();
		list.add("Extractor");
		list.add("a=" + in + "Arabidopsis_lyrata.v.1.0.31.chr.gff3");
		list.add("g=" + in + "Arabidopsis_lyrata.v.1.0.31.dna.genome.fa");
		list.add("r=" + true );
		list.add("f=" + false );
		list.add("Ambiguity=" + Ambiguity.AMBIGUOUS );
		list.add("p=" + true );
					
		String[] given = {
				"cds-parts.fasta",
				"assignment.tabular",
				"proteins.fasta",
		};
		
		cliCheck(list, given, given);
	}
	
	@Test
	public void testERE() throws Exception {
		ArrayList<String> list = new ArrayList<String>();
		list.add("ERE");
		list.add("s="+Stranded.FR_UNSTRANDED );
		list.add("m=" + in + "RNAseq.bam");
		list.add("c="+ true );

		String[] given = {
				"coverage.bedgraph",
				"introns.gff"
		};
		
		cliCheck(list, given, given);
	}
	
	@Test
	public void testGeMoMa() throws Exception {
		ArrayList<String> list = new ArrayList<String>();
		list.add("GeMoMa");
		
		boolean simple = true;
		//only first part
		//list.add("selected=" + in + "selected-first.txt" );
		
		list.add("t=" + in + "tblastn.tabular");
		list.add("tg=" + in + "TAIR10_chr_all.fas");
		list.add("a=" + in + "assignment.tabular");
		list.add("c=" + in + "cds-parts.fasta");
		list.add("q=" + in + "proteins.fasta");
		if( simple ) {
			list.add("i=" + in + "introns.gff");
		} else {
			list.add("i=" + in + "at-introns-1.gff");
			list.add("i=" + in + "at-introns-2.gff");
		}
		if( simple ) {
			list.add("coverage=UNSTRANDED");
			list.add("coverage_unstranded=" + in + "coverage.bedgraph");
		} else {
			list.add("coverage=UNSTRANDED");
			list.add("coverage_unstranded=" + in + "coverage.bedgraph-0.txt");
			list.add("coverage=UNSTRANDED");
			list.add("coverage_unstranded=" + in + "coverage.bedgraph-1.txt");
		}
		
		String[] given = {
				"predicted_annotation.gff",
				"predicted_protein.fasta"
		};
		
		cliCheck(list, given, given);
	}
	
	@Test
	public void testGAF() throws Exception {
		ArrayList<String> list = new ArrayList<String>();
		list.add("GAF");
		
		list.add("g=" +in+ "predicted_annotation.gff");
		
		String[] given = {"filtered_predictions.gff"};
		
		cliCheck(list, given, given);
	}
}
