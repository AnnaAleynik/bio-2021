require 'pry'

DNA_COMPLEMENTARY = {
  "a" => "t",
  "t" => "a",
  "c" => "g",
  "g" => "c"
}

GENETIC_TABLE = {
  "TTT" => "Phe", "TCT" => "Ser", "TAT" => "Tyr",  "TGT" => "Cys",
  "TTC" => "Phe", "TCC" => "Ser", "TAC" => "Tyr",  "TGC" => "Cys",
  "TTA" => "Leu", "TCA" => "Ser", "TAA" => "STOP", "TGA" => "STOP",
  "TTG" => "Leu", "TCG" => "Ser", "TAG" => "STOP", "TGG" => "Trp",
  
  "CTT" => "Leu", "CCT" => "Pro", "CAT" => "His", "CGT" => "Arg",
  "CTC" => "Leu", "CCC" => "Pro", "CAC" => "His", "CGC" => "Arg",
  "CTA" => "Leu", "CCA" => "Pro", "CAA" => "Gln", "CGA" => "Arg",
  "CTG" => "Leu", "CCG" => "Pro", "CAG" => "Gln", "CGG" => "Arg",
  
  "ATT" => "Ile", "ACT" => "Thr", "AAT" => "Asn", "AGT" => "Ser",
  "ATC" => "Ile", "ACC" => "Thr", "AAC" => "Asn", "AGC" => "Ser",
  "ATA" => "Ile", "ACA" => "Thr", "AAA" => "Lys", "AGA" => "Arg",
  "ATG" => "Met", "ACG" => "Thr", "AAG" => "Lys", "AGG" => "Arg",
  
  "GTT" => "Val", "GCT" => "Ala", "GAT" => "Asp", "GGT" => "Gly",
  "GTC" => "Val", "GCC" => "Ala", "GAC" => "Asp", "GGC" => "Gly",
  "GTA" => "Val", "GCA" => "Ala", "GAA" => "Glu", "GGA" => "Gly",
  "GTG" => "Val", "GCG" => "Ala", "GAG" => "Glu", "GGG" => "Gly"
}.freeze

NO_FOUND_MESSAGE = "No orf"
DNA_LENGTH_FORMAT_ERROR = "Enter length of DNA from 100 to 1000."
GC_CONTENT_FORMAT_ERROR = "Enter GC content from 20 to 80."

STOP_CODONS = ["tga", "taa", "tag"].freeze
START_CODONS = ["atg"].freeze

def input
  return @input if @input

  puts "Enter length of DNA:"
  dna_length = gets.chomp.to_i
  # if NOT number was entered it will be 0
  abort(DNA_LENGTH_FORMAT_ERROR) if dna_length < 100 || dna_length > 1000

  puts "Enter GC content:"
  gc_content_procent = gets.chomp.to_i
  # 1/2 - C, 1/2 - G, = A + T
  # GC=60 => G – 30%, C – 30%, A – 20%, T – 20%
  abort(GC_CONTENT_FORMAT_ERROR) if gc_content_procent < 20 || gc_content_procent > 80

  @input ||= { dna_length: dna_length, gc_content_procent: gc_content_procent }
end

def find_orf(dna)
  max_length = 0
  start_found = false
  start_index = -1
  curr_start_index = -1
  # protein

  # print "frame "
  codons = dna.each_slice(3).map(&:join)
  codons.each_with_index do |val, i|
    if START_CODONS.include?(val) && !start_found
      # binding.pry
      # print "start #{i}: #{val} "
      curr_start_index = i
      start_found = true
    end

    if STOP_CODONS.include?(val) && start_found
      l = i - curr_start_index + 1
      # max_length = max_length < l ? l : max_length
      if max_length < l
        max_length = l
        start_index = curr_start_index
      end
      start_found = false
      # print "length: #{l} i:#{i}"
    end
  end

  # print " st #{start_index}"
  # puts "max_length #{max_length}"
  max_length = max_length * 3

  { length: max_length, start: start_index}
end

def explore_three_frames(dna, direct)
  frames = [find_orf(dna)]

  first_letter = dna.shift
  frames << find_orf(dna)

  second_letter = dna.shift
  frames << find_orf(dna)

  dna.unshift(first_letter, second_letter)

  frames.each_with_index do |f, i|
    f[:start] = f[:start] * 3 + (i + 1)
    f[:frame] = i + 1
    f[:direct] = direct
  end
end

def find_max_orf(frames)
  # frames.each { |f| puts f}

  frames.delete_if { |frame| frame[:length] < 30 }

  abort(NO_FOUND_MESSAGE) if frames.empty?

  frames.reduce do |f, o|
    f[:length] > o[:length] ? f : o
  end
end

def orf_seq(orf, dna, reversed_dna)
  if orf[:direct] == "direct"
    # puts " in orf_seg #{orf[:start]}"
    dna[orf[:start] - 1, orf[:length]]
  else
    # puts " in orf_seg #{orf[:start]}"
    reversed_dna[orf[:start] - 1, orf[:length]]
  end
end

# begin

dna_length = input[:dna_length]
gc_content_procent = input[:gc_content_procent]
# dna_length = 150
# gc_content_procent = 45

g_number = dna_length * gc_content_procent / 100 / 2
c_number = g_number
a_number = (dna_length - g_number - c_number) / 2
t_number = a_number

dna = ("a" * a_number + "t" * t_number + "g" * g_number + "c" * c_number).chars.shuffle
# dna = "ttcgccatgacgggtctgatgtgggactcgcctatctggtactcatgggtgaacgccgaggttatctacctctaagcgcgccgaattcgcaacggcaggtagtacccacgcctcataggctccgtacccctgagaaccctttcagttaggtaacggaagcgcggtagtgtcggtcccagcgggggtgacaacatgccgttacaaccattgtttcgtagtacggtgtggcacatccacaggttcgccgcgctggagtcagtagacgggctggtagcctcgcctcggatcggtcttactagtatgcacagctccacaaagcgctccctactcttgtcggcgcgtatcgttggtcggacaggaagacgcattgggggtgcacgagcccagatatggtagacttccggagcctacgcgtccgatttagcacagtaacgcggtgggaccgagccggctccccagaccacgccgtcccaactccgattacggcaagtagcaccagctctcctaccggaacccatcgcatggttcatgtgtgggttaaggtgctatgcagcacatcgtactaggggctcaatctccgaacgggtatccgcatgctcaggggatacgagccgcgctcgtcgcgggggtatgggcccggacggcaccctaagcactaaaacgcggctcttcccagattagtagagcacgcctaggcgccttcgacccgggccccgctccggatcagacaggacggggagcctgatgtcgagtcgcatggcacccccgtcgtcggtgtgaacggagattaatgtgtcgcc".chars
puts "DNA: #{dna.join}"

reversed_dna = dna.join.gsub(Regexp.union(DNA_COMPLEMENTARY.keys), DNA_COMPLEMENTARY).reverse.chars

found_frame = [
  explore_three_frames(dna, "direct"),
  explore_three_frames(reversed_dna, "reverse")
].flatten!

orf = find_max_orf(found_frame)
orf_seq = orf_seq(orf, dna, reversed_dna)

protein = orf_seq.each_slice(3).map(&:join).map { |codon| GENETIC_TABLE[codon.upcase] }.join(" ")

puts "ORF: #{orf_seq.each_slice(3).map(&:join).join(" ")}"
# protein
puts "Protein: #{protein}"
# range
puts "From #{orf[:start]} to #{orf[:start] + orf[:length] - 1}"
# reverse or direct
puts "ORF is found in the #{orf[:direct]} strand."
# frame
puts "ORF is found in the #{orf[:frame]} frame."
