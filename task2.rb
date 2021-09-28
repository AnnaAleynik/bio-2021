require 'zlib'
require 'benchmark'

puts "Enter file name:"
current_path = File.dirname(__FILE__)
file_name = gets.chomp
file_n = File.new(current_path + "/test_data/" + file_name)

i = 0
count = 0
average_length = 0
quality_percent = 0


# Benchmark.bm(1) do |x|
#   x.report do
    # Zlib::GzipReader.open("./test_data/test_example3.fastq.gz") do |file|
    Zlib::GzipReader.open(file_n) do |file|
      file.each do |entry|
        print "." if count % 5000 == 0

        count += 1
        i += 1

        if i == 2
          average_length += entry.length - 1
        elsif i == 4
          entry.each_char do |c|
            quality_percent += 1 if c < "J" && c > ">"
          end
          i = 0
        end
      end
    end
#   end
# end
puts

reads_count = count / 4
letters_count = average_length

puts "count         #{reads_count}"
puts "letters count #{letters_count}"
puts "average       #{(average_length.to_f/reads_count).round(1)}"
puts "quality       #{(quality_percent.to_f * 100/letters_count).round(2)}%"
