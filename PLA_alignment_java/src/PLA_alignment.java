/**
 * Java program for PLA alignment from raw reads
 * Command line argument example: input="abc\def.txt" (no whitespaces are allowed around the equal signs)
 */


// Import packages
import java.io.*;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.TreeSet;
import java.util.Objects;
import java.util.zip.*;
import java.util.Random;

//import java.lang.Math;
import java.time.LocalDate;
import java.time.LocalTime;
import java.time.LocalDateTime;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.text.similarity.LevenshteinDistance;
import com.google.common.collect.*;


public class PLA_alignment
{
	
	// Calculate Hamming distance
	// Can allow N, which matches to all bases
	public static int HammingDistanceCalculator(String str1, String str2, boolean allowN)
	{
		if (str1.length() != str2.length()) {throw new java.lang.IllegalArgumentException("Input strings must have equal lengths!");}
		else
		{
			int dist = 0;
			
			for (int i=0; i<str1.length(); i++)
			{
				
				if (str1.charAt(i) != str2.charAt(i) )
				{
					if (allowN && ((str1.charAt(i)=='N') || (str2.charAt(i)=='N'))) {continue;} else {dist++;}
				}
			}
			return dist;
		}
	}
	public static int HammingDistanceCalculator(String str1, String str2) {return HammingDistanceCalculator(str1, str2, false);} // overloading: default allowN is false
	
	public static void main(String[] args)
	{
		
		switch (args[0])
		{
		
		/**
		 * Alignment from raw reads for drop-seq runs
		 * R1=... R2=... O=... ABfile=... SUMMARY=... HEADER=... (for doing alignment of paired reads 1 by 1)
		 * 		or
		 * dirI=... dirO=... SUMMARY=... (for doing alignment of all paired reads contained in directory dirI)
		 * Input arguments:
		 * 		R1: path to read 1 file (fastq.gz format)
		 * 		R2: path to read 2 file (fastq.gz format)
		 * 		ABfile: path to PLA target-DNA barcode lookup table (csv format)
		 * 		O: path to store output file (txt.gz format)
		 * 		SUMMARY: directory to store summary files (txt format) (default is current working directory)
		 * 		HEADER: whether the ABfile has header to be skipped (default is false)
		 * 	or
		 * 		dirI: path to the directory containing read 1 and read 2 files (fastq.gz format)
		 * 		dirO: path to the directory to store output file (txt.gz) ---> the output file names are taken from the files in dirI
		 * 		SUMMARY: directory to store summary files (txt format)
		 * 
		 * Output format: cell barcodes , UMI , AB1 ID , AB2 ID
		 */
		case "ReadAlignmentDropSeq":
		{
			
			// Parse the arguments
			String R1 = "", R2 = "", ABfile = "", O = "";
			String SUMMARY = System.getProperty("user.dir") + File.separator + "ReadAlignmentDropSeq_summary.txt"; // default summary file directory
			boolean skip_header=false;
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "R1": R1 = j[1]; break;
				case "R2": R2 = j[1]; break;
				case "ABfile": ABfile = j[1]; break;
				case "O": O = j[1]; break;
				case "SUMMARY": SUMMARY = j[1]; break;
				case "HEADER": skip_header = "true".equalsIgnoreCase(j[1]); break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
		
			try (
					BufferedReader br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(R1)))); // read1
					BufferedReader br2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(R2)))); // read2
					BufferedReader brAB = new BufferedReader(new FileReader(ABfile)); // AB-DNA barcode look up table
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))); // output file
					BufferedWriter bwsum = new BufferedWriter(new FileWriter(SUMMARY)) // summary file
				 )
			{
				
				// Write out the command line arguments
				System.out.println();
				for (String j : args)
				{
					bwsum.write(j); bwsum.newLine();
					System.out.println(j);
				}
				bwsum.newLine();
				System.out.println();
				
				// Read the AB look up table into an array
				List<List<String>> ABarray = new ArrayList<List<String>>();
				String ABline;
				if (skip_header) {brAB.readLine();} // skip first line
				while ((ABline=brAB.readLine()) != null)
				{
					String[] values = ABline.split(",");
					ABarray.add(Arrays.asList(values));
				}
	
				
				/**
				 * Process read 1 and 2
				 * Read 1 contains cell barcode and UMI
				 * Read 2 contains PLA products
				 * Output format: <cell barcode> , <UMI> , <AB1 ID> , <AB2 ID>
				 */
				
				// Set up the counters for the summary files
				int short_read_counter = 0; // counter for the number of reads with fewer than 75 bases
				int badUMI_counter = 0; // counter for the number of reads with invalid UMIs due to repeated G
				int excessiveG_counter = 0; // counter for the number of reads with too many Ns
				int excessiveN_counter = 0; // counter for the number of reads with too many Ns
				int bad_connector_counter = 0; // counter for the number of reads with non-matching connector sequence
				int non_matching_AB_counter = 0; // counter of the number of reads with non-matching AB barcode
				
				// Set up for alignment
				String line1, line2;
				int counter = 0; // line number counter
				int PLA_counter = 0; // PLA product counter
				String connector = "TCGTGTCGTGTCGTGTCTAAAG"; // connector sequence
				
				// Start reading
				long my_timer = System.currentTimeMillis();
				System.out.printf("%10s   %-8s   ReadAlignmentDropSeq   Start alignment%n", LocalDate.now(), LocalTime.now().withNano(0));
				bwsum.write("Start alignment at " + LocalDateTime.now()); bwsum.newLine();
				
				while (((line1=br1.readLine()) != null) && ((line2=br2.readLine()) != null))
				{
							
					if ((counter % 4) == 1)
					{
						
						
						
						// Expected location of connector starts at index 39 (base 40th)
						int connector_start = 39;					
						int connector_start_temp = 39; // temporary starting location of connector, for used in for loop
						
						
						// Skip reads with at least 7 occurrences of G in the UMI region
						if (line2.length() < 75)
						{
							short_read_counter++;
						}
						
						else if (StringUtils.countMatches(line1.substring(12), "G") >= 7)
						{
							badUMI_counter++;
						}
						
						// Skip reads with excessive Gs (reads that contain only G in the connector region)
						else if ( ( StringUtils.countMatches(line2.substring(connector_start,connector_start+connector.length()), "G") ) == connector.length() )
						{
							excessiveG_counter++;
						}
						
						// Check if read 2 has more than 10 N's
						else if (StringUtils.countMatches(line2, "N") > 10)
						{
							excessiveN_counter++;
						}
						
						else
						{							
							// Locate the connector using Levenshtein distance, max allowed distance is 3
							// Does not need to take into account N, since the connector region doesn't usually contain N
							int[] connector_shift = new int[] {0, -1, 1};
							boolean match_connector = false; // true if matching connector is found
							
							int connector_distance = 2; // lowest Levenshtein distance found between the true connector and the test connector sequence
							for (int shift_i : connector_shift)
							{
								int temp_distance = LevenshteinDistance.getDefaultInstance().apply(line2.substring(connector_start+shift_i, connector_start+shift_i+connector.length()), connector);
								if ((temp_distance <= 2) && (temp_distance < connector_distance))
								{
									connector_distance = temp_distance;
									connector_start_temp = connector_start + shift_i;
									match_connector = true;
									
									if (temp_distance == 0)
									{ break; }
								}
							}
							connector_start = connector_start_temp;
							
							if (match_connector)
							{
								
								// Initialize the AB ID
								String AB1_ID = "Unknown";
								String AB2_ID = "Unknown";
							
								// Initialize variables to store hamming distance for each AB1 and 2 barcode
								int[] temp1 = new int[ABarray.size()];
								int[] temp2 = new int[ABarray.size()];
								
								// Check if there is a frameshift, in order to locate AB2 correctly
								int shift_j = line2.substring(connector_start, connector_start+connector.length()+1).indexOf("TAAAG"); // location of AAAG in the found connector region
								if (shift_j == -1)
								{
									counter++;
									bad_connector_counter++;
									continue;
								}
								shift_j = 17 - shift_j; // number of bases to shift to the left
								
								// Found AB barcodes
								String AB1_found = line2.substring(connector_start-18, connector_start-18+8);
								String AB2_found = line2.substring(connector_start+25-shift_j, connector_start+25-shift_j+8);
								
								// Calculate Hamming distance only if there is at most 1 N in the found barcode
								if ((StringUtils.countMatches(AB1_found, "N") <= 1) && (StringUtils.countMatches(AB2_found, "N") <= 1))
								{
									// Counter for # of matches with 1 hamming distance
									int match_counter1 = 0, match_counter2 = 0;
									// Index for matches with 1 hamming distance
									int match_index1 = -1, match_index2 = -1;
									
									for (int i=0; (i<ABarray.size()) && ((AB1_ID=="Unknown") || (AB2_ID=="Unknown")) && ((match_counter1<2) && (match_counter2<2)); i++)
									{
										// Calculate Hamming distance
										temp1[i] = HammingDistanceCalculator(AB1_found, ABarray.get(i).get(1), true);
										temp2[i] = HammingDistanceCalculator(AB2_found, ABarray.get(i).get(1), true);
										
										// Allow early termination of for loop if found an exact match
										if (temp1[i] == 0)
										{
											AB1_ID = ABarray.get(i).get(0);
										}
										else if (temp1[i] == 1) {match_counter1++; match_index1 = i;}
										
										if (temp2[i] == 0)
										{
											AB2_ID = ABarray.get(i).get(0);
										}
										else if (temp2[i] == 1) {match_counter2++; match_index2 = i;}
										
									}
									
									// Find unambiguous match with 1 Hamming distance (ie, discard reads that have more than 1 matches with 1 hamming distance
									if (match_counter1 == 1) {AB1_ID = ABarray.get(match_index1).get(0);}
									if (match_counter2 == 1) {AB2_ID = ABarray.get(match_index2).get(0);}
									
									if (!Objects.equals(AB1_ID,"Unknown") && !Objects.equals(AB2_ID,"Unknown"))
									{
										bwout.write(line1.substring(0, 12)+","+line1.substring(12)+","+AB1_ID+","+AB2_ID);
										bwout.newLine();
										PLA_counter++;
									}
									else
									{
										non_matching_AB_counter++;
									}
								}
							}
							
							else
							{
								bad_connector_counter++;
							}
							
							
							
						}
						
						
						
						
						if ((((counter-1)/4+1) % 1000000) == 0)
						{
							System.out.printf("%10s   %-8s   ReadAlignmentDropSeq   Processed %,15d records   Elapsed time for last 1,000,000 records: %ds%n",
									LocalDate.now(), LocalTime.now().withNano(0), (counter-1)/4+1, (System.currentTimeMillis()-my_timer)/1000);
							my_timer = System.currentTimeMillis();
						}
					}
					
					
//					if ((((counter-1)/4)+1)>2000000) {System.out.printf("%s   %s   Processed %,d lines%n", LocalDate.now(), LocalTime.now().withNano(0), counter); break;} // for testing purposes
				
					counter++;
				}
				
				
				System.out.printf("%10s   %-8s   ReadAlignmentDropSeq   Done: processed %,d records%n", LocalDate.now(), LocalTime.now().withNano(0), (counter-1)/4+1);
				System.out.printf("\tNumber of valid PLA products: %,15d%n", PLA_counter);
				
				// Write to summary file
				bwsum.write("ReadAlignmentDropseq: Finished at " + LocalDateTime.now().withNano(0) + ", processed " + String.format("%,d",(counter-1)/4+1) + " records"); bwsum.newLine();
				bwsum.write("Number of valid PLA products: " + String.format("%,d", PLA_counter)); bwsum.newLine();
				bwsum.write("Number of records discarded because of read 2 being too short: " + String.format("%,d",short_read_counter)); bwsum.newLine();
				bwsum.write("Number of records discarded because of excessive G in UMIs: " + String.format("%,d",badUMI_counter)); bwsum.newLine();
				bwsum.write("Number of records discarded because of excessive number of Ns: " + String.format("%,d",excessiveN_counter)); bwsum.newLine();
				bwsum.write("Number of records discarded because of excessive G in read 2: " + String.format("%,d",excessiveG_counter)); bwsum.newLine();
				bwsum.write("Number of records discarded because of non-matching connector sequence: " + String.format("%,d",bad_connector_counter)); bwsum.newLine();
				bwsum.write("Number of records discarded because of non-matching antibody barcode: " + String.format("%,d",non_matching_AB_counter)); bwsum.newLine();

			} catch (IOException e) { throw new IllegalArgumentException("Invalid file paths!");}
			break;
			}
			
		/**
		 * Alignment from raw reads for smart-seq runs
		 * R1=... O=... ABfile=... SUMMARY=... HEADER=... DOWNSAMPLE=... SEED=...
		 * 
		 * Input arguments:
		 * 		R1: path to read 1 file (fastq.gz format)
		 * 		O: path to store output file (txt.gz format)
		 * 		ABfile: path to PLA target-DNA barcode lookup table (csv format)
		 * 		SUMMARY: directory to store summary files (txt format) (default is current working directory)
		 * 		HEADER: whether the ABfile has header to be skipped (default is false)
		 * 		DOWNSAMPLE: downsample parameter (default is 1)
		 * 		^^^^^^^^^^ if 0 < DOWNSAMPLE < 1, it is the downsample ratio (1 = no downsample, 0.6 = use 60% of the reads)
		 * 		^^^^^^^^^^ if 1 < DOWNSAMPLE, it is the number of reads to retain (10,000 = use the first 10,000 reads)
		 * 		SEED: integer seed number for random downsampling or bootstrap (default = 1)
		 * 
		 * Output format: cell barcodes, UMI , AB1 ID , AB2 ID
		 * 
		 * The summary file also contains the found AB barcode and their read counts (only keep those with read counts >= 100)
		 * NOTE: actual seed number for random downsampling is actually SEED + DOWNSAMPLE*1000, so that different downsample ratios don't use the same seed number even though the SEED values are equal
		 */
		case "ReadAlignmentSmartSeq":
		{
			
			// Parse the arguments
			String R1 = "", ABfile = "", O = "";
			double DOWNSAMPLE = 1;
			int SEED = 1;
			String SUMMARY = System.getProperty("user.dir") + File.separator + "ReadAlignmentSmartSeq_summary.txt"; // default summary file directory
			boolean skip_header=false;
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "R1": R1 = j[1]; break;
				case "ABfile": ABfile = j[1]; break;
				case "O": O = j[1]; break;
				case "SUMMARY": SUMMARY = j[1]; break;
				case "HEADER": skip_header = "true".equalsIgnoreCase(j[1]); break;
				case "DOWNSAMPLE": DOWNSAMPLE = Double.parseDouble(j[1]); break;
				case "SEED": SEED = Integer.parseInt(j[1]); break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
		
			try (
					BufferedReader br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(R1)))); // read1
					BufferedReader brAB = new BufferedReader(new FileReader(ABfile)); // AB-DNA barcode look up table
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))); // output file
					BufferedWriter bwsum = new BufferedWriter(new FileWriter(SUMMARY)) // summary file
				 )
			{
				
				// Write out the command line arguments
				System.out.println();
				for (String j : args)
				{
					bwsum.write(j); bwsum.newLine();
					System.out.println(j);
				}
				bwsum.newLine();
				System.out.println();
			
				// Get the cell barcode from the file name
				String[] filename = R1.split("/");
				filename = filename[filename.length-1].split("\\.");
				String cellBC = filename[0];
				
				// Read the AB look up table into an array
				List<List<String>> ABarray = new ArrayList<List<String>>();
				String ABline;
				if (skip_header) {brAB.readLine();} // skip first line
				while ((ABline=brAB.readLine()) != null)
				{
					String[] values = ABline.split(",");
					ABarray.add(Arrays.asList(values));
				}
	
				// Hash Multiset to store unique UMIs and their number of occurrences
				Multiset<String> ABcounts = HashMultiset.create();
				
				/**
				 * Process read 1
				 * Read 1 contains PLA products
				 * UMI region: 2nd base to 17th base (16-base long)
				 * Output format: <cell barcode> , <UMI> , <AB1 ID> , <AB2 ID>
				 * Cell barcode is equal to sample barcode
				 */

				// Set up the counters for the summary files
				int counter = 0; // line number counter
				int downsample_counter = 0; // counters of processed reads (for downsampling)
				int PLA_counter = 0; // PLA product counter
				int bad_connector_counter = 0; // counter for reads with non-matching connector
				int bad_UMI_counter = 0; // counter of reads with bad UMI region
				int non_matching_AB_counter = 0; // counter of the number of reads with non-matching AB barcode
				
				// Set up for alignment
				String line1;
				String connector = "TCGTGTCGTGTCGTGTCTAAAG"; // connector sequence
				Random generator = new Random(SEED + (int) DOWNSAMPLE*1000); // Initialize a random number generator
				
				// If using downsample with a specific number of reads, get the number of reads first
//				List<Integer> downsample_read = new ArrayList<Integer>();
//				int temp_counter = 0;
//				while ((line1=br1.readLine()) != null)
//				{
//					downsample_read.add(temp_counter);
//					temp_counter++;
//				}
//				br1.close();
//				br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(R1)))); // read1
				
				// Start reading
				long my_timer = System.currentTimeMillis();
				System.out.printf("%10s   %-8s   Start alignment%n", LocalDate.now(), LocalTime.now().withNano(0));
				bwsum.write("Start alignment at " + LocalDateTime.now()); bwsum.newLine();
				
				while ((line1=br1.readLine()) != null)
				{
							
					if ((counter % 4) == 1)
					{				
						
						// Downsample if required
						if ((DOWNSAMPLE < 1) && (DOWNSAMPLE > 0))
						{
							if (generator.nextDouble() > DOWNSAMPLE)
							{
								if ((((counter-1)/4+1) % 1_000_000) == 0)
								{
									System.out.printf("%10s   %-8s   ReadAlignmentSmartSeq   Processed %,15d records   Elapsed time for last 1,000,000 records: %ds%n",
											LocalDate.now(), LocalTime.now().withNano(0), (counter-1)/4+1, (System.currentTimeMillis()-my_timer)/1000);
									my_timer = System.currentTimeMillis();
								}
								
								counter++;
								continue;
							}
						}
						else if (DOWNSAMPLE > 1)
						{
							if ((counter-1)/4 >= DOWNSAMPLE)
							{
								break;
							}
						}
						else if (DOWNSAMPLE < 0)
						{
							throw new java.lang.IllegalArgumentException("Invalid downsample parameter!");
						}
						
						// Increment the counter for the number of processed reads
						downsample_counter++;
						
						// Check if the UMI region has an excessive amount of A
						if (StringUtils.countMatches(line1.substring(1, 17), "A") >= 10)
						{
							if ((((counter-1)/4+1) % 1_000_000) == 0)
							{
								System.out.printf("%10s   %-8s   ReadAlignmentSmartSeq   Processed %,15d records   Elapsed time for last 1,000,000 records: %ds%n",
										LocalDate.now(), LocalTime.now().withNano(0), (counter-1)/4+1, (System.currentTimeMillis()-my_timer)/1000);
								my_timer = System.currentTimeMillis();
							}
							
							counter++;
							bad_UMI_counter++;
							continue;
						}
						
						int connector_start = 39; // expected starting location of connector is at index 39 (base 40th)
						int connector_start_temp = 39; // temporary starting location of connector, for used in for loop
						
						// Locate the connector using Levenshtein distance, max allowed distance is 2
						// Does not need to take into account N, since the connector region doesn't usually contain N
						int[] connector_shift = new int[] {0, -1, 1, -2, 2, -3, 3}; // only allow up to 3-base shift of the expected starting location
						boolean match_connector = false; // true if matching connector is found
						int connector_distance = 3; // lowest Levenshtein distance found between the true connector and the test connector sequence
						for (int shift_i : connector_shift)
						{
							int temp_distance = LevenshteinDistance.getDefaultInstance().apply(line1.substring(connector_start+shift_i, connector_start+shift_i+connector.length()), connector);
							if ((temp_distance <= 2) && (temp_distance < connector_distance))
							{
								connector_distance = temp_distance;
								connector_start_temp = connector_start + shift_i;
								match_connector = true;
								
								if (temp_distance == 0)
								{ break; }
							}
						}
						connector_start = connector_start_temp;
						
						if (match_connector)
						{
							
							// Initialize the AB ID
							String AB1_ID = "Unknown";
							String AB2_ID = "Unknown";
						
							// Initialize variables to store hamming distance for each AB1 and 2 barcode
							int[] temp1 = new int[ABarray.size()];
							int[] temp2 = new int[ABarray.size()];
							
							// Check if there is a frameshift, in order to locate AB2 correctly
							// Find the location of TAAAG in the found connector region
							int shift_j = line1.substring(connector_start, connector_start+connector.length()+2).indexOf("TAAAG"); // add 2 just in case there are 2 insertions
							if (shift_j == -1) // skip read if can't find TAAAG in the connector region
							{
								counter++;
								bad_connector_counter++;
								continue;
							}
							shift_j = 17 - shift_j; // number of bases to shift to the left
							
							// Found AB barcodes
							String AB1_found = line1.substring(connector_start-18, connector_start-18+8);
							String AB2_found = "";
							if ((connector_start+25-shift_j+8) <= line1.length()) // check if the read fully contains the AB2 ID
							{
								AB2_found = line1.substring(connector_start+25-shift_j, connector_start+25-shift_j+8);
							}
							else if ((connector_start+25-shift_j+8) == (line1.length()+1)) // the first 7 bases of the barcode is at the end of the read (in other words, some insertions)
							{
								AB2_found = line1.substring(connector_start+25-shift_j) + "N";
							}
							else // skip read if there are too many insertions
							{
								counter++;
								continue;
							}
							
							
							
							// Calculate Hamming distance only if there is at most 1 N in the found barcode
							if ((StringUtils.countMatches(AB1_found, "N") <= 1) && (StringUtils.countMatches(AB2_found, "N") <= 1))
							{
								// Counter for # of matches with 1 hamming distance
								int match_counter1 = 0, match_counter2 = 0;
								// Index for matches with 1 hamming distance
								int match_index1 = -1, match_index2 = -1;

								
								for (int i=0; (i<ABarray.size()) && (Objects.equals(AB1_ID,"Unknown") || Objects.equals(AB2_ID,"Unknown")) && ((match_counter1<=1) && (match_counter2<=1)); i++)
								{
									// Calculate Hamming distance
									temp1[i] = HammingDistanceCalculator(AB1_found, ABarray.get(i).get(1), true);
									temp2[i] = HammingDistanceCalculator(AB2_found, ABarray.get(i).get(1), true);
									
									// Allow early termination of for loop if found an exact match
									if (Objects.equals(AB1_ID,"Unknown"))
									{
										if (temp1[i] == 0)
										{
											AB1_ID = ABarray.get(i).get(0);
										}
										else if (temp1[i] == 1)
										{match_counter1++; match_index1 = i;}
									}
									
									if (Objects.equals(AB2_ID,"Unknown"))
									{
										if (temp2[i] == 0)
										{
											AB2_ID = ABarray.get(i).get(0);
										}
										else if (temp2[i] == 1)
										{match_counter2++; match_index2 = i;}
									}
									
								}
								
								// Find unambiguous match with 1 Hamming distance (ie, discard reads that have more than 1 matches with 1 hamming distance
								if (match_counter1 == 1) {AB1_ID = ABarray.get(match_index1).get(0);}
								if (match_counter2 == 1) {AB2_ID = ABarray.get(match_index2).get(0);}
								
								// There can be a deletion between the connector and the end of AB1 ID --> Check for this
								if ((match_counter1 == 0) && Objects.equals(AB1_ID,"Unknown"))
								{
									AB1_found = line1.substring(connector_start-18+1, connector_start-18+1+8);
									match_counter1 = 0;
									match_index1 = -1;
									
									for (int i=0; (i<ABarray.size()) && Objects.equals(AB1_ID,"Unknown") && (match_counter1<=1); i++)
									{
										// Calculate Hamming distance
										temp1[i] = HammingDistanceCalculator(AB1_found, ABarray.get(i).get(1), true);
										if (temp1[i] == 0)
										{
											AB1_ID = ABarray.get(i).get(0);
										}
										else if (temp1[i] == 1)
										{match_counter1++; match_index1 = i;}
									}
									
									if (match_counter1 == 1) {AB1_ID = ABarray.get(match_index1).get(0);}
								}
								
								if (!Objects.equals(AB1_ID,"Unknown") && !Objects.equals(AB2_ID,"Unknown"))
								{
									bwout.write(cellBC + "," + line1.substring(1, 17) + "," + AB1_ID + "," + AB2_ID);
									bwout.newLine();
									PLA_counter++;
								}
								else
								{
									non_matching_AB_counter++;
								}
							}
							
							// Add the found AB barcodes to the Hash Multiset ABcounts
							ABcounts.add(AB1_found);
							ABcounts.add(AB2_found);
							
						}
						else
						{
							bad_connector_counter++;
						}
						
						
						
						if ((((counter-1)/4+1) % 1_000_000) == 0)
						{
							System.out.printf("%10s   %-8s   ReadAlignmentSmartSeq   Processed %,15d records   Elapsed time for last 1,000,000 records: %ds%n",
									LocalDate.now(), LocalTime.now().withNano(0), (counter-1)/4+1, (System.currentTimeMillis()-my_timer)/1000);
							my_timer = System.currentTimeMillis();
						}
					}
					
					
//					if ((((counter-1)/4)+1)>2000000) {System.out.printf("%s   %s   Processed %,d lines%n", LocalDate.now(), LocalTime.now().withNano(0), counter); break;} // for testing purposes
				
					counter++;
				}
				
				
				System.out.printf("%10s   %-8s   ReadAlignmentSmartSeq   Done: processed %,d reads%n", LocalDate.now(), LocalTime.now().withNano(0), (counter+1)/4);
				System.out.printf("\tNumber of reads with a valid PLA product: %,15d%n", PLA_counter);
				
				// Write to summary file
				bwsum.write("ReadAlignmentSmartSeq: Finished at " + LocalDateTime.now().withNano(0) + ", processed " + String.format("%,d",(counter-1)/4+1) + " reads"); bwsum.newLine();
				bwsum.write("Number of reads that are processed: " + String.format("%,d", downsample_counter)); bwsum.newLine();
				bwsum.write("Number of reads with a valid PLA product: " + String.format("%,d", PLA_counter)); bwsum.newLine();
				bwsum.write("Number of reads discarded because of non-matching connector sequence: " + String.format("%,d",bad_connector_counter)); bwsum.newLine();
				bwsum.write("Number of reads discarded because of bad UMI: " + String.format("%,d",bad_UMI_counter)); bwsum.newLine();
				bwsum.write("Number of records discarded because of non-matching antibody barcode: " + String.format("%,d",non_matching_AB_counter)); bwsum.newLine();

				// Add to the summary file the found AB barcodes
				bwsum.newLine();
				bwsum.write("Antibody barcode\tOccurrences"); bwsum.newLine();
				// Sort the HashMultiset by decreasing occurrences, and save the top 10 to an array
				String[] AB_sortedbycounts = Multisets.copyHighestCountFirst(ABcounts).elementSet().toArray(new String[0]);
				int temp_counter = 0;
				for (String i : AB_sortedbycounts)
				{
					// Break for loop once read count is lower than 100
					if (temp_counter >= 10)
					{break;}
					
					bwsum.write(i + String.format("\t%,10d", ABcounts.count(i)));
					bwsum.newLine();
					temp_counter++;
				}
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths!");}
			break;
			}
		

		/**
		 * Cell barcode correction from aligned reads (ie, output of PLAReadAlignment)
		 * 		I=... O=... SUMMARY=...
		 * 
		 * Input arguments:
		 * 		I: path to aligned reads (txt.gz format)
		 * 		O: path to store output file (txt.gz format)
		 * 		SUMMARY: directory to store summary files (txt format) (default is current working directory)
		 * 		CELL_BC_LIST: path to a comma-separated list of cell barcodes produced by drop-seq pipeline (rows are cell barcodes, column 0 is readcount, column 1 is the cell barcode sequence)
		 * 		^^^^^^^^^^^^ output of drop-seq tools' BAMTagHistogram function (txt.gz format)
		 * 		READCOUNT_CUTOFF: only keep the barcode sequence with at least this number of readcount (default is 1000)
		 * 		HEADER: whether the CELL_BC_LIST have header, which will be skipped (default is false)
		 * 
		 * Correction method: n-gram with Hamming distance <=1
		 * 		Split the query cell barcode into 2 halves at the first N
		 * 		Example: barcode "123N567N9" is split into "123" and "567N9"
		 * 		If there is no N, check for n-grams split in the middle: "123456" is split into "123" and "456"
		 * 		Then check for matching "123" and "5678" in the reference n-grams
		 * 		Look for perfectly matched cell barcodes, or cell barcodes that are unambiguously mismatched by 1 Hamming distance (ie, discard cell barcodes that have more than 1 matching cell barcodes with Hamming = 1) 
		 */	
		case "CellBarcodeCorrection":
		{
			
			// Parse the arguments
			String I = "", O = "", cell_BC_path = "";
			String SUMMARY = System.getProperty("user.dir") + File.separator + "CellBarcodeCorrection_summary.txt"; // default summary file directory
			int rc_cutoff = 1000;
			boolean skip_header = false;
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "I": I = j[1]; break;
				case "O": O = j[1]; break;
				case "SUMMARY": SUMMARY = j[1]; break;
				case "CELL_BC_LIST": cell_BC_path = j[1]; break;
				case "READCOUNT_CUTOFF": rc_cutoff = Integer.parseInt(j[1]); break;
				case "HEADER": skip_header = "true".equalsIgnoreCase(j[1]); break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			try (
					BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(I)))); // aligned reads
					BufferedReader brBC = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(cell_BC_path)))); // cell barcode list
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))); // output file
					BufferedWriter bwsum = new BufferedWriter(new FileWriter(SUMMARY)) // summary file
				)
			{
				
				// Write out the command line arguments
				for (String j : args)
				{
					bwsum.write(j + " "); bwsum.newLine();
					System.out.println(j);
				}
				bwsum.newLine();
				System.out.println();
				
				
				// Set up the counters for the summary files
				int PLA_counter = 0; // counter for the number of valid PLA reads
				int exact_counter = 0; // counter for the number of reads with exact match
				int ambiguous_counter = 0; // counter for the number of reads that matches ambiguously with reference cell barcodes (ie, more than 1 cell barcode with 1 hamming distance)
				int nomatch_counter = 0; // counter for the number of reads with no match
				int manyN_counter = 0; // counter for the number of reads with more than 1 N in the cell barcode region
				
				// Skip first line if required
				if (skip_header) {brBC.readLine();}
				
//				// Initialize a multimap of key-value pair: key is the sequence, and value is an ArrayList of the sequence's index in the reference barcode list
				ListMultimap<String, Integer> cellBCmap_1 = ArrayListMultimap.create(); // first segment of the barcode
				ListMultimap<String, Integer> cellBCmap_2 = ArrayListMultimap.create(); // second segment of the barcode
				
				// Read accepted cell barcode sequence into an ArrayList
				System.out.printf("%10s   %-8s   CellBarcodeCorrection   Start processing reference cell barcodest%n", LocalDate.now(), LocalTime.now().withNano(0));
				List<String> cellBCarray = new ArrayList<String>(); // ArrayList containing the cell barcode
				List<String> cellBCarray_1 = new ArrayList<String>(); // ArrayList containing the first half of the cell barcode (base 1-6)
				List<String> cellBCarray_2 = new ArrayList<String>(); // ArrayList containing the second half of the cell barcode (base 7-12)
				
				// Set up counters for reference cell barcodes
				int ref_counter = 0; // counter for the number of cell barcodes
				
				// Read the cell barcodes and build the multimap for the non-N case
				String BCline;
				while ((BCline=brBC.readLine()) != null)
				{
					String[] values = BCline.split("\t");
					if (Integer.parseInt(values[0]) >= rc_cutoff)
					{
						cellBCarray.add(values[1]);
						// Split the reference barcode into 2 equal halves, and store them in the ArrayList
						cellBCarray_1.add(values[1].substring(0, values[1].length()/2));
						cellBCarray_2.add(values[1].substring(values[1].length()/2));
						
						// Store the 2 halves in the multimap
						cellBCmap_1.put(cellBCarray_1.get(cellBCarray_1.size()-1), ref_counter);
						cellBCmap_2.put(cellBCarray_2.get(cellBCarray_2.size()-1), ref_counter);
						ref_counter++;
					}
				}
				
				// Split the reference barcode at different positions for the N case
				// 		Example: barcode "abcdef", N can be at any position
				// Initialize a list of multimaps of key-value pair
				// Each multimap is the n-gram and its index in the original reference list, split at a specific location
				//		Example: the 2nd multimap in the list contains the n-gram split at the 2nd base: "aNcdef"
				// If N is at the beginning or the end of the cell barcode, ignore N, and split at the middle: "Nbcdef" is split into "bc" and "def", "abcdeN" is split into "abc" and "de"
				List<ListMultimap<String, Integer>> cellBCmaplist_1 = new ArrayList<>();
				List<ListMultimap<String, Integer>> cellBCmaplist_2 = new ArrayList<>();
				
				// Initialize an ArrayList of ArrayList
				// Each sub-ArrayList is all the segments of the reference cell barcode (non-collasped), split at a specific location
				List<ArrayList<String>> cellBCarraylist_1 = new ArrayList<>();
				List<ArrayList<String>> cellBCarraylist_2 = new ArrayList<>();
				
				System.out.printf("%10s   %-8s   CellBarcodeCorrection   Start building maps of reference cell barcodes%n", LocalDate.now(), LocalTime.now().withNano(0));
				int barcode_length = cellBCarray.get(0).length();
				for (int i=0; i<barcode_length; i++) // loop through each possible splitting position
				{
					// Temporary variables
					ListMultimap<String, Integer> tempmap_1 = ArrayListMultimap.create();
					ListMultimap<String, Integer> tempmap_2 = ArrayListMultimap.create();
					ArrayList<String> temparray_1 = new ArrayList<>();
					ArrayList<String> temparray_2 = new ArrayList<>();
					
					// Specify start index (inclusive) of the first n-gram, and end index (exclusive) of the second n-gram
					int temp_start = 0;
					int temp_end = barcode_length;
					// Specify split position
					int temp_split_1 = i;
					int temp_split_2 = i + 1;
					if (i == 0)
					{
						temp_start = 1;
						temp_split_1 = barcode_length/2;
						temp_split_2 = barcode_length/2;
					}
					else if (i == (barcode_length-1))
					{
						temp_end = barcode_length - 1;
						temp_split_1 = barcode_length/2;
						temp_split_2 = barcode_length/2;
					}
					
					for (int j=0; j<cellBCarray.size(); j++) // loop through reference cell barcodes
					{
						String temp_1 = cellBCarray.get(j).substring(temp_start, temp_split_1);
						String temp_2 = cellBCarray.get(j).substring(temp_split_2, temp_end);
						tempmap_1.put(temp_1, j);
						tempmap_2.put(temp_2, j);
						
						temparray_1.add(temp_1);
						temparray_2.add(temp_2);
					}
					
					cellBCmaplist_1.add(tempmap_1);
					cellBCmaplist_2.add(tempmap_2);
					cellBCarraylist_1.add(temparray_1);
					cellBCarraylist_2.add(temparray_2);
				}
				
				System.out.printf("%10s   %-8s   CellBarcodeCorrection   Finish reading %,d cell barcodes%n", LocalDate.now(), LocalTime.now().withNano(0), ref_counter);
				
				
				// n-gram method: takes into account N
				String line;
				int counter = 0;
				long my_timer = System.currentTimeMillis();
				System.out.printf("%10s   %-8s   Start cell barcode correction%n", LocalDate.now(), LocalTime.now().withNano(0));
				bwsum.write("Start cell barcode correction at " + LocalDateTime.now()); bwsum.newLine();
				while ((line=brI.readLine()) != null)
				{
						
					// Time stamp (put here because there are continue statements below)
					if ((counter>0) && (counter % 1000000 == 0))
					{
						System.out.printf("%10s   %-8s   CellBarcodeCorrection   Processed %,15d records   Elapsed time for last 1,000,000 records: %ds%n",
								LocalDate.now(), LocalTime.now().withNano(0), counter, (System.currentTimeMillis()-my_timer)/1000);
						my_timer = System.currentTimeMillis();
					}
					
					// Read the cell barcodes, UMI and AB IDs
					List<String> values = new ArrayList<String>(Arrays.asList(line.split(",")));
					
					// counter for how many barcodes have hamming distance <= 1, discard cell barcode if counter >=2
					int match_counter = 0;
					
					// Case 1: query barcode doesn't contain N
					if (values.get(0).indexOf("N") == -1)
					{
						// Split cell barcodes into half, and compare each to the reference cell barcodes
						String testBC_1 = values.get(0).substring(0, values.get(0).length()/2);
						String testBC_2 = values.get(0).substring(values.get(0).length()/2);
						
						// Check if the first half has a match
						if (cellBCmap_1.containsKey(testBC_1))
						{
							// Look for matching second half
							for (int i : cellBCmap_1.get(testBC_1))
							{
								if (Objects.equals(testBC_2, cellBCarray_2.get(i)))
								{
									exact_counter++;
									match_counter = 1;
									break;
								}
								else if (HammingDistanceCalculator(testBC_2,cellBCarray_2.get(i)) == 1)
								{
									match_counter++;
									values.set(0, cellBCarray.get(i));
								}
							}
						}
						// Check if the second half has a match
						else if (cellBCmap_2.containsKey(testBC_2))
						{
							// Look for matching first half with 1 hamming distance
							for (int i : cellBCmap_2.get(testBC_2))
							{
								if (HammingDistanceCalculator(testBC_1, cellBCarray_1.get(i)) == 1)
								{
									match_counter++;
									values.set(0, cellBCarray_1.get(i) + testBC_2);
								}
								if (match_counter >= 2)
								{
									break;
								}
								
							}
						}
					}
					
					// Case 2: query barcode does contain N
					else
					{
						// Allow up to 1 N's only
						if (StringUtils.countMatches(values.get(0), "N") > 1)
						{
							manyN_counter++;
							counter++;
							continue;
						}
						
						int split_position = values.get(0).indexOf("N"); // position of N
						
						// Split query cell barcode into n-gram
						String testBC_1, testBC_2;
						if (split_position == 0)
						{
							testBC_1 = values.get(0).substring(1, values.get(0).length()/2);
							testBC_2 = values.get(0).substring(values.get(0).length()/2 );
						}
						else if (split_position == (values.get(0).length()-1))
						{
							testBC_1 = values.get(0).substring(0, values.get(0).length()/2);
							testBC_2 = values.get(0).substring(values.get(0).length()/2, values.get(0).length() - 1);
						}
						else
						{
							testBC_1 = values.get(0).substring(0, split_position);
							testBC_2 = values.get(0).substring(split_position + 1);
						}
						
						
						// Check if the first half has a match
						if (cellBCmaplist_1.get(split_position).containsKey(testBC_1))
						{
							for (int i : cellBCmaplist_1.get(split_position).get(testBC_1))
							{
								if (Objects.equals(testBC_2, cellBCarraylist_2.get(split_position).get(i)))
								{
									exact_counter++;
									match_counter = 1;
									values.set(0, cellBCarray.get(i));
									break;
								}
								else if (HammingDistanceCalculator(testBC_2, cellBCarraylist_2.get(split_position).get(i)) == 1)
								{
									match_counter++;
									values.set(0, cellBCarray.get(i));
								}
							}
						}
						
						else if (cellBCmaplist_2.get(split_position).containsKey(testBC_2))
						{
							for (int i : cellBCmaplist_2.get(split_position).get(testBC_2))
							{
								if (HammingDistanceCalculator(testBC_1, cellBCarraylist_1.get(split_position).get(i)) == 1)
								{
									match_counter++;
									values.set(0, cellBCarray.get(i));
								}
								
								if (match_counter >= 2)
								{
									break;
								}
							}
						}
					}
					
								
					if (match_counter == 1)
					{
						PLA_counter++;
						bwout.write(String.join(",", values));
						bwout.newLine();
					}
					else if (match_counter == 0) // no matching barcode is found
					{
						nomatch_counter++;
					}
					else if (match_counter >= 2)
					{
						ambiguous_counter++;
					}
					
					counter++;
					
				}
				
				System.out.printf("%10s   %-8s   CellBarcodeCorrection   Done%n", LocalDate.now(), LocalTime.now().withNano(0));
				System.out.printf("\tNumber of valid PLA products: %,15d%n", PLA_counter);
				
				// Write to summary file
				bwsum.write("CellBarcodeCorrection: Finished at " + LocalDateTime.now().withNano(0) + ", processed " + String.format("%,d",counter) + " records"); bwsum.newLine();
				bwsum.write("Number of accepted PLA products: " + String.format("%,d", PLA_counter)); bwsum.newLine();
				bwsum.write("Number of exact matches: " + String.format("%,d",exact_counter)); bwsum.newLine();
				bwsum.write("Number of discarded ambiguous reads: " + String.format("%,d",ambiguous_counter)); bwsum.newLine();
				bwsum.write("Number of discarded non-matching reads: " + String.format("%,d",nomatch_counter)); bwsum.newLine();
				bwsum.write("Number of discarded cell barcodes with more than 1 N base: " + String.format("%,d", manyN_counter)); bwsum.newLine();
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths!");}
			break;
		}

		
		/**
		 * Check if any UMIs map to more than 1 PLA product
		 * 		I=... O=...
		 * Input arguments:
		 * 		I: path to aligned/cell barcode-corrected file (txt.gz format)
		 * 		O: path to output multi-mapping UMIs (tab-separated, txt.gz format)
		 */
		
		case "CheckUMIMapping":
		{
			// Parse the arguments
			String I = "", O = "";
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "I": I = j[1]; break;
				case "O": O = j[1]; break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			try (
					BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(I)))); // aligned reads
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))) // output file
				  )
			{
				// Write out the command line arguments
				for (String j : args)
				{
					System.out.println(j);
				}
				System.out.printf("%n");
				
				// Initialize a SetMultimap
				// keys: unique cell barcode-UMI combinations (sorted alphabetically)
				// values: set of corresponding unique PLA products
				SetMultimap<String, String> UMI_multimap = MultimapBuilder.treeKeys().hashSetValues().build();
				
				// Time stamp
				long my_timer = System.currentTimeMillis();
				System.out.printf("%10s   %-8s   Start CheckUMIMapping%n", LocalDate.now(), LocalTime.now().withNano(0));
				
				// Add the aligned reads to the multimap
				String line;
				while ((line = brI.readLine()) != null)
				{
					String[] values = line.split(",");
					UMI_multimap.put(values[0] + "_" + values[1], values[2] + ":" + values[3]);
				}
				
				// Export duplicated UMI mappings
				int counter = 0;
				Set<String> temp_values = new HashSet<>();
				for (String i : UMI_multimap.keySet())
				{
					temp_values = UMI_multimap.get(i);
					if (temp_values.size() > 1)
					{
						bwout.write(i + "\t" + temp_values);
						bwout.newLine();
					}
					counter++;
					
					if ((counter % 100_000) == 0)
					{
						System.out.printf("%10s   %-8s   CheckUMIMapping   Processed %,15d UMIs   Elapsed time for last 100,000 UMIs: %ds%n",
								LocalDate.now(), LocalTime.now().withNano(0), (counter-1)/4+1, (System.currentTimeMillis()-my_timer)/1000);
						my_timer = System.currentTimeMillis();
					}
				}
				
				// Time stamp
				System.out.printf("%10s   %-8s   CheckUMIMapping   Done: processed %,d UMIs%n", LocalDate.now(), LocalTime.now().withNano(0), counter);
				
			} catch (IOException e) { throw new IllegalArgumentException("Invalid file paths!");}
			
			break;
		}
		
		
		/**
		 * Perform UMI merging on cell barcode-corrected PLA products
		 * 		I=... O=... SUMMARY=...
		 * Input arguments:
		 * 		I: path to aligned/cell barcode-corrected file (txt.gz format)
		 * 		O: path to output UMI merged file (txt.gz format)
		 * 		SUMMARY: path to store summary file
		 * 
		 * Merge the following cases:
		 * 		Exact matches
		 * 		Match with 1 hamming distance to another UMI that has a higher read count
		 * 		(This means that barcodes that match at 1 hamming distance to more than 1 other barcodes are discarded)
		 */
		
		case "UMIMerging":
		{
			// Parse the arguments
			String I = "", O = "";
			String SUMMARY = System.getProperty("user.dir") + File.separator + "UMIMerging_summary.txt"; // default summary file directory
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "I": I = j[1]; break;
				case "O": O = j[1]; break;
				case "SUMMARY": SUMMARY = j[1]; break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			try (
					BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(I)))); // aligned reads
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))); // output file
					BufferedWriter bwsum = new BufferedWriter(new FileWriter(SUMMARY)) // summary file
				  )
			{

				// Write out the command line arguments
				for (String j : args)
				{
					bwsum.write(j); bwsum.newLine();
					System.out.println(j);
				}
				bwsum.newLine();
				System.out.printf("%n");
				
				// Each PLA product is described by 2 features: [cell barcode/AB1 ID/AB2 ID], and [UMI]. We shall call these feature 1 and 2, respectively
				// Initialize
				List<String> UMIarray = new ArrayList<String>(); // ArrayList to store all UMI
				ListMultimap<String, Integer> PLAproductmap = ArrayListMultimap.create(); // ArrayListMultimap to store unique feature 1's and its index (for indexing UMIarray)
				
				// Read the file
				System.out.printf("%10s   %-8s   UMIMerging   Start pre-processing records%n", LocalDate.now(), LocalTime.now().withNano(0));
				bwsum.write("Start UMIMerging at " + LocalDateTime.now().withNano(0)); bwsum.newLine();
				long my_timer = System.currentTimeMillis();
				String line;
				int record_counter = 0;
				while ((line = brI.readLine()) != null)
				{
					String[] values = line.split(",");
					UMIarray.add(values[1]);
					PLAproductmap.put(values[0] + "," + values[2] + "," + values[3], record_counter);
					record_counter++;
					if ((record_counter % 1_000_000) == 0)
					{
						System.out.printf("%10s   %-8s   UMIMerging   Processed %,15d records   Elapsed time for last 1,000,000 records: %ds%n",
								LocalDate.now(), LocalTime.now().withNano(0), record_counter, (System.currentTimeMillis()-my_timer)/1000);
						my_timer = System.currentTimeMillis();
					}
				}
				
				System.out.printf("%10s   %-8s   UMIMerging   Finished pre-processing%n", LocalDate.now(), LocalTime.now().withNano(0));
				System.out.printf("\tNumber of unique cell barcode-PLA pair combinations: %,15d%n", PLAproductmap.keySet().size());
				
				// Loop through each unique feature 1, and merge UMI from the same feature 1 that are within 1 hamming distance
				int counter = 0; // total record counter
				int unique_counter = 0; // unique PLA product counter
				my_timer = System.currentTimeMillis();
				System.out.printf("%10s   %-8s   UMIMerging   Start merging UMI%n", LocalDate.now(), LocalTime.now().withNano(0));
				for (String str_i : PLAproductmap.keySet())
				{

					// Hash Multiset to store unique UMIs and their number of occurences
					Multiset<String> temp_UMIcounts = HashMultiset.create();
					// Hashset to store unique UMIs, from which duplicated UMIs will be removed
					Set<String> temp_UMI = new HashSet<>();
					
					// Get the UMIs associated with this feature 1
					for (int j : PLAproductmap.get(str_i))
					{
						temp_UMIcounts.add(UMIarray.get(j));
						temp_UMI.add(UMIarray.get(j));
					}
					
					// Hashset to store unique UMIs, from which duplicated UMIs have been removed
					Set<String> temp_UMI_updated = new HashSet<>(temp_UMI);
					
					// Sort the HashMultiset by decreasing occurrences, and save to an array
					String[] temp_UMI_sortedbycounts = Multisets.copyHighestCountFirst(temp_UMIcounts).elementSet().toArray(new String[0]);
					
					// Iterate through the unique UMIs in order of increasing occurrences, and remove UMIs that are within 1 hamming distance with at least one other, more frequent, UMI
					// Slow step**** //
					for (int j=(temp_UMI_sortedbycounts.length-1); j>=0; j--)
					{
						temp_UMI = new HashSet<>(temp_UMI_updated);
						for (String str_j : temp_UMI)
						{
							if (HammingDistanceCalculator(str_j, temp_UMI_sortedbycounts[j], true) == 1)
							{
								temp_UMI_updated.remove(temp_UMI_sortedbycounts[j]);
								break;
							}
						}
					}
					
					// Write out the UMI corrected reads
					String[] values = str_i.split(",");
					for (String str_j : temp_UMI_updated)
					{
						bwout.write(values[0] + "," + str_j + "," + values[1] + "," + values[2]); bwout.newLine();
						unique_counter++;
					}
					
					// Time stamp
					counter++;
					if ((counter % 1_000) == 0)
					{
						System.out.printf("%10s   %-8s   UMIMerging   Processed %,12d records   Elapsed time for last 1,000 records: %ds%n",
								LocalDate.now(), LocalTime.now().withNano(0), counter, (System.currentTimeMillis()-my_timer)/1000);
						my_timer = System.currentTimeMillis();
					}
				}
				
				System.out.printf("%10s   %-8s   UMIMerging   Done%n", LocalDate.now(), LocalTime.now().withNano(0));
				System.out.printf("\tNumber of unique PLA products: %,15d%n", unique_counter);
				
				// Write to summary file
				bwsum.write("UMIMerging: Finished at " + LocalDateTime.now().withNano(0) + ", processed " + String.format("%,d",counter) + " unique cell barcode-PLA pair combinations"); bwsum.newLine();
				bwsum.write("Number of unique PLA products: " + String.format("%,d", unique_counter)); bwsum.newLine();

				
				
			} catch (IOException e) { throw new IllegalArgumentException("Invalid file paths!");}
			break;
		}
		
		
		/**
		 * Tally the read counts each cell barcode receives
		 * 		I=... O=...
		 * Input arguments
		 * 		I: path to cell barcode corrected file (txt.gz format)
		 * 		O: path to store output file (txt.gz format)
		 * Merge the following cases:
		 * 		Exact matches
		 * 		Match with 1 hamming distance to exactly one other UMI
		 * This means that barcodes that match at 1 hamming distance to more than 1 other sequences are discarded 
		 */
		
		case "ReadcountHistogram":
		{
			// Parse the arguments
			String I = "", O = "";
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "I": I = j[1]; break;
				case "O": O = j[1]; break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			try (
					BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(I)))); // aligned reads
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))) // output file
				)
			{
				
				System.out.printf("%10s   %-8s   ReadcountHistogram   Start%n", LocalDate.now(), LocalTime.now().withNano(0));
				
				// Create Hash Multiset to get counts of each cell barcode
				Multiset<String> readcounts = HashMultiset.create();
				
				// Read file
				String line;
				while ( (line = brI.readLine()) != null )
				{
					String[] values = line.split(",");
					readcounts.add(values[0]);
				}
				
				// Write to output file, tab-separated
				// First row is header, first column is read count, second column is cell barcode
				bwout.write("Readcounts\tCellbarcode"); bwout.newLine();
				for (String i : Multisets.copyHighestCountFirst(readcounts).elementSet())
				{
					bwout.write(String.format("%d", readcounts.count(i)) + "\t" + i); bwout.newLine();
				}
				
				System.out.printf("%10s   %-8s   ReadcountHistogram   Done%n", LocalDate.now(), LocalTime.now().withNano(0));
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths!");}
			
			break;
		}
		
		
		/**
		 * Get digital count from UMI merged file, output a tab-separated txt.gz file
		 * 		I=... O=... CELL_BC_LIST=... SUMMARY=...
		 * Input arguments:
		 * 		I: path to the UMI merged file (txt.gz format)
		 * 		O: path to store the digital count matrix (txt.gz format)
		 * 		CELL_BC_LIST: a list of chosen cell barcodes (from knee plot) (txt format)
		 * 		HEADER: whether the CELL_BC_LIST has header to be skipeed (default is false)
		 * 		SUMMARY: path to store the summary file
		 */
		
		case "DigitalCountDropSeq":
		{

			
			// Parse the arguments
			String I = "", O = "", cell_BC_path = "";
			boolean skip_header = false;
			String SUMMARY = System.getProperty("user.dir") + File.separator + "DigitalCountDropSeq_summary.txt"; // default summary file directory
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "I": I = j[1]; break;
				case "O": O = j[1]; break;
				case "CELL_BC_LIST": cell_BC_path = j[1]; break;
				case "HEADER": skip_header = "true".equalsIgnoreCase(j[1]); break;
				case "SUMMARY": SUMMARY = j[1]; break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			try (
					BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(I)))); // aligned reads
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))); // output file
					BufferedReader brBC = new BufferedReader(new FileReader(cell_BC_path)); // cell barcode list
					BufferedWriter bwsum = new BufferedWriter(new FileWriter(SUMMARY)) // summary file
				)
			{
				// Write out the command line arguments
				for (String j : args)
				{
					bwsum.write(j); bwsum.newLine();
					System.out.println(j);
				}
				bwsum.newLine();
				System.out.printf("%n");
				
				// Time stamp
				System.out.printf("%10s   %-8s   DigitalCountDropSeq   Start%n", LocalDate.now(), LocalTime.now().withNano(0));
				
				// Read the list of chosen cell barcodes
				// hash set to store chosen cell barcodes, ordered by read counts from drop-seq tools (0 index = most reads)
				Set<String> chosen_BC = new LinkedHashSet<>();
				String BCline;
				if (skip_header) {brBC.readLine();}
				while ((BCline = brBC.readLine()) != null)
				{
					chosen_BC.add(BCline);
				}
				
				// Read the PLA products, only keep chosen cell barcodes
				String line;
				Multiset<String> PLAproduct = HashMultiset.create(); // count the UMIs for each combination of cell barcode-PLA pair
				Set<String> PLA_ID = new TreeSet<>(); // hash set to store unique PLA pairs, ordered by AB1 ID alphabetically
				while ((line = brI.readLine()) != null)
				{
					String[] values = line.split(",");
					if (chosen_BC.contains(values[0]))
					{
						PLAproduct.add(values[0] + "_" + values[2] + ":" + values[3]);
						PLA_ID.add(values[2] + ":" + values[3]);
					}
				}
				
				// Write the digital count matrix row by row (ie, by PLA pair)
				// Header: PLA_pair \t cell barcode 1 \t cell barcode 2 ...
				bwout.write("PLA_pair");
				for (String str_i : chosen_BC)
				{
					bwout.write("\t" + str_i);
				}
				bwout.newLine();
				
				// Write each row
				for (String str_i : PLA_ID)
				{
					bwout.write(str_i);
					for (String str_j : chosen_BC)
					{
						bwout.write(String.format("\t%d", PLAproduct.count(str_j + "_" + str_i)));
					}
					bwout.newLine();
				}
				
				// Time stamp
				System.out.printf("%10s   %-8s   DigitalCountSmartSeq   Done%n", LocalDate.now(), LocalTime.now().withNano(0));
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths!");}
						
			break;
		}
		
		/**
		 * Get digital count from UMI merged files for Smart-Seq data, output a tab-separated txt.gz file
		 * 		dirI=... O=... ABfile=... SUMMARY=... HEADER=... REMOVE_DUPLICATE=... EXPORT_DUPLICATE=...
		 * Input arguments:
		 * 		dirI: path to the directory containing ONLY the UMI merged files (txt.gz format)
		 * 		O: path to store the digital count matrix (txt.gz format) ***** Rows are samples, columns are PLA pairs
		 * 		ABfile: path to PLA target-DNA barcode lookup table (csv format)
		 * 		SUMMARY: path to store the summary file
		 * 		HEADER: whether the ABfile has header to be skipped (default is false)
		 * 		REMOVE_DUPLICATE: whether to remove duplicated PLA products across cells (ie, PLA molecules that have the same UMI and PLA IDs) (default is false)
		 *		EXPORT_DUPLICATE: whether to export the full list of duplicated PLA products (txt.gz format) (default is false)
		 *		^^^^^^^^^^^^^^^^ the export file is saved in the current working directory as export_duplicate.txt.gz
		 */
		
		case "DigitalCountSmartSeq":
		{

			// Parse the arguments
			String dirI = "", O = "", ABfile = "";
			boolean skip_header = false;
			String SUMMARY = System.getProperty("user.dir") + File.separator + "DigitalCountSmartSeq_summary.txt"; // default summary file directory
			boolean remove_duplicate = false;
			boolean export_duplicate = false;
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "dirI": dirI = j[1]; break;
				case "O": O = j[1]; break;
				case "ABfile": ABfile = j[1]; break;
				case "SUMMARY": SUMMARY = j[1]; break;
				case "HEADER": skip_header = "true".equalsIgnoreCase(j[1]); break;
				case "REMOVE_DUPLICATE": remove_duplicate = "true".equalsIgnoreCase(j[1]); break;
				case "EXPORT_DUPLICATE": export_duplicate = "true".equalsIgnoreCase(j[1]); break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			try (
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))); // output file
					BufferedReader brAB = new BufferedReader(new FileReader(ABfile)); // AB-DNA barcode look up table
					BufferedWriter bwsum = new BufferedWriter(new FileWriter(SUMMARY)) // summary file
				)
			{
				// Write out the command line arguments
				for (String j : args)
				{
					bwsum.write(j); bwsum.newLine();
					System.out.println(j);
				}
				bwsum.newLine();
				System.out.printf("%n");
				
				System.out.printf("%10s   %-8s   DigitalCountSmartSeq   Start%n", LocalDate.now(), LocalTime.now().withNano(0));
				
				// Read the list of possible AB IDs
				List<String> AB_list = new ArrayList<>();
				String ABline;
				if (skip_header) {brAB.readLine();}
				while ((ABline = brAB.readLine()) != null)
				{
					AB_list.add(ABline.split(",")[0]);
				}
				
				// Make an ArrayList of all possible PLA pairs
				List<String> PLA_ID = new ArrayList<>();
				for (String AB1_ID : AB_list)
				{
					for (String AB2_ID : AB_list)
					{
						PLA_ID.add(AB1_ID + ":" + AB2_ID);
					}
				}
				
				// Write the first row: header
				bwout.write("sample_name");
				for (String str_i : PLA_ID)
				{
					bwout.write("\t" + str_i);
				}
				bwout.newLine();
				
				// Prepare the variables for reading the input directory
				File dir = new File(dirI);
				File[] dir_list = dir.listFiles();
				
				// If REMOVE_DUPLICATE = true, build a list of duplicated PLA products
				Multiset<String> PLAduplicate = HashMultiset.create(); // count the occurrences of each unique UMI+PLA product
				if ((remove_duplicate) || (export_duplicate))
				{
					System.out.printf("%10s   %-8s   DigitalCountSmartSeq   Building the list of duplicated PLA products%n", LocalDate.now(), LocalTime.now().withNano(0));
					if (dir_list != null)
					{
						for (File child : dir_list)
						{
							// Skip non .txt.gz files
							if (!child.getName().endsWith(".txt.gz"))
							{
								continue;
							}
							
							try (BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(child)))))
							{
								// Read the PLA products and store the counts of each PLA product across all cells in the HashMultiset
								String line;
								while ((line = brI.readLine()) != null)
								{
									// Split the reads into 2 part: cell barcode + PLA product
									String[] values = line.split(",", 2);
									
									// Add to the HashMultiset PLAduplicate
									PLAduplicate.add(values[1]);
								}
							} catch (IOException e) {throw new IllegalArgumentException("Invalid file format in the input directory!");}
						}
					}
				}
				
				// Read the files in directory dirI
				int counter = 0; // counter for number of processed files
				int PLA_counter = 0; // counter for number of total PLA products (before removing duplicates)
				int duplicate_counter = 0; // counter for the number of counts that have been removed
				long my_timer = System.currentTimeMillis();
				
				System.out.printf("%10s   %-8s   DigitalCountSmartSeq   Writing the DGE file%n", LocalDate.now(), LocalTime.now().withNano(0));
				if (dir_list != null)
				{
					for (File child : dir_list)
					{
						try (BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(child)))))
						{
							
							// Read the PLA products and get the UMI counts
							String line;
							String sample_name = "";
							Multiset<String> PLAproduct = HashMultiset.create(); // count the UMIs for each PLA pair (PLA ID)
							while ((line = brI.readLine()) != null)
							{
								String[] values = line.split(",");
								sample_name = values[0];
								PLA_counter++;
								
								if (remove_duplicate && (PLAduplicate.count(values[1] + "," + values[2] + "," + values[3]) > 1)) // check for duplicated PLA product
								{
									duplicate_counter++;
								}
								else
								{
									PLAproduct.add(values[2] + ":" + values[3]);
								}
								
							}
							
							// Write the digital count matrix
							// Header: sample_name \t PLA_ID 1 \t PLA_ID 2 ...
							bwout.write(sample_name);
							for (String str_i : PLA_ID)
							{
								bwout.write(String.format("\t%d", PLAproduct.count(str_i)));
							}
							bwout.newLine();
							
							
						} catch (IOException e) {throw new IllegalArgumentException("Invalid file format in the input directory!");}
						
						
						counter++;
						if ((counter % 20) == 0)
						{
							System.out.printf("%10s   %-8s   ReadAlignment   Processed %,15d files   Elapsed time for last 20 files: %ds%n",
									LocalDate.now(), LocalTime.now().withNano(0), counter, (System.currentTimeMillis()-my_timer)/1000);
							my_timer = System.currentTimeMillis();
						}
					}
				}
				else
				{
					throw new IllegalArgumentException("Empty file directory!");
				}
				
				// Time stamp
				System.out.printf("%10s   %-8s   DigitalCountSmartSeq   Done%n", LocalDate.now(), LocalTime.now().withNano(0));
				
				// Record the number of removed duplicated counts
				bwsum.write("Number of total PLA products across all samples: " + String.format("%,d",PLA_counter)); bwsum.newLine();
				bwsum.write("Number of duplicated counts that were removed across all samples: " + String.format("%,d",duplicate_counter));	bwsum.newLine();
				
				// Record top 10 duplicated PLA products
				bwsum.write("PLA product\tNumber of occurences"); bwsum.newLine();
				int temp_counter = 0;
				for (String i : Multisets.copyHighestCountFirst(PLAduplicate).elementSet())
				{
					if (temp_counter >= 10)
					{
						break;
					}
					bwsum.write(String.format("%s\t%d", i, PLAduplicate.count(i)));
					bwsum.newLine();
						
					temp_counter++;
				}
				
				// Export list of all PLA duplicates
				if (export_duplicate)
				{
					String export_duplicate_file = System.getProperty("user.dir") + File.separator + "export_duplicate.txt.gz";
					try (BufferedWriter bwdup = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(export_duplicate_file)))))
					{
						bwdup.write("PLA product\tNumber of occurences"); bwdup.newLine();
						for (String i : Multisets.copyHighestCountFirst(PLAduplicate).elementSet())
						{
							if (PLAduplicate.count(i) == 1)
							{
								break;
							}
							bwdup.write(String.format("%s\t%d", i, PLAduplicate.count(i))); bwdup.newLine();
						}
					} catch (IOException e) {;}
				}
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths!");}
						
			break;
		}
		
		default: throw new java.lang.IllegalArgumentException("Invalid function call!");
		}
	}

}