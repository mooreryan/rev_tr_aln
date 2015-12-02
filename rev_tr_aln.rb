#!/usr/bin/env ruby

# Copyright 2015 Ryan Moore
# Contact: moorer@udel.edu
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

Signal.trap("PIPE", "EXIT")

require "set"
require "parse_fasta"
require "trollop"

module Utils
  def self.time fmt="%F %T.%L"
    Time.now.strftime fmt
  end

  def self.log msg, *args, &block
    $stderr.puts "[INFO] -- " +
                 "#{msg % args}"

    block.call if block
  end

  def self.warn msg, *args, &block
    $stderr.puts "[WARN] -- " +
                 "#{msg % args}"

    block.call if block
  end

    def self.assert test, msg, *args, &block
    unless test
      err_msg = "[ERROR] (#{self.caller}) -- #{msg}"

      if block
        block.call
      end

      abort err_msg % args
    end
  end

  def self.check_file fname
    self.assert File.exists?(fname),
                "File '%s' does not exist.",
                fname
  end

  def self.caller start=2
    Kernel.caller(start).join ", "
  end

  def self.parse_fname fname
    { dir: File.dirname(fname),
      base: File.basename(fname, File.extname(fname)),
      ext: File.extname(fname) }
  end

  def self.trollop_check_file arg, name
    if arg.nil?
      Trollop.die name, "You didn't provide an input file!"
    elsif !File.exists?(arg)
      Trollop.die name, "#{arg} doesn't exist!"
    end

    self.parse_fname arg
  end
end

def get_clean_hash fname
  orfs = {}
  FastaFile.open(fname).each_record do |head, seq|
    orfs[head] = seq.upcase
  end

  orfs
end

def get_acc header
  header.split(" ").first
end

def degap seq
  seq.gsub(/[-\.]/, "")
end

def char_len dna
  dna.gsub(/[-\.]/, "").length
end

DNA_TO_AA = {}

AA_TO_DNA = {
  "I" => %w[ATT ATC ATA],
  "L" => %w[CTT CTC CTA CTG TTA TTG],
  "V" => %w[GTT GTC GTA GTG],
  "F" => %w[TTT TTC],
  "M" => %w[ATG],
  "C" => %w[TGT TGC],
  "A" => %w[GCT GCC GCA GCG],
  "G" => %w[GGT GGC GGA GGG],
  "P" => %w[CCT CCC CCA CCG],
  "T" => %w[ACT ACC ACA ACG],
  "S" => %w[TCT TCC TCA TCG AGT AGC],
  "Y" => %w[TAT TAC],
  "W" => %w[TGG],
  "Q" => %w[CAA CAG],
  "N" => %w[AAT AAC],
  "H" => %w[CAT CAC],
  "E" => %w[GAA GAG],
  "D" => %w[GAT GAC],
  "K" => %w[AAA AAG],
  "R" => %w[CGT CGC CGA CGG AGA AGG],
  "*" => %w[TAA TAG TGA]
}

AA_SET = Set.new AA_TO_DNA.keys

opts = Trollop.options do
  banner <<-EOS

  The first bit is the acc number. Headers in protein alignment and
  prot orf file must be the same.

  Codon table from http://www.cbs.dtu.dk/courses/27619/codon.html

  Options:
  EOS

  opt(:protein_alignment, "Protein alignment file", type: :string)
  opt(:aa_orfs, "AA orfs", type: :string)
  # opt(:translated_nucleotides, "6 frame translation of nt orfs",
  #     type: :string)
  opt(:nt_orfs, "NT orfs", type: :string)
  # opt(:outdir, "Output directory", type: :string, default: ".")
end

Utils.trollop_check_file(opts[:protein_alignment], :protein_alignment)
Utils.trollop_check_file(opts[:aa_orfs], :aa_orfs)
Utils.trollop_check_file(opts[:nt_orfs], :nt_orfs)

# aa_head => [nt_head1, nt_head2, ... ]
ident_seqs = {}
tr_matches = 0
uniq_tr_matches = 0

nt_orfs = get_clean_hash opts[:nt_orfs]
aa_orfs = get_clean_hash opts[:aa_orfs]
# trans_orfs = get_clean_hash opts[:translated_nucleotides]

# check for similar headers
aa_orfs.each do |aa_head, aa_seq|
  nt_orfs.each do |nt_head, nt_seq|
    aa_acc = get_acc aa_head
    nt_acc = get_acc nt_head

    if aa_acc == nt_acc
      if ident_seqs.has_key? aa_head
        # TODO can this ever actually happen?
        ident_seqs[aa_head] << nt_head
      else
        ident_seqs[aa_head] = [nt_head]
      end
    end
  end
end

Utils.log "Found %d acc matches", ident_seqs.count

Utils.log "%d AA seqs have multiple matches",
          ident_seqs.select { |k, v| v.count > 1 }.count

still_need_nt_seqs = aa_orfs.keys - ident_seqs.keys

# still_need_nt_seqs.each do |aa_head|
#   aa_seq = aa_orfs[aa_head]
#   Utils.assert aa_seq, "%s not found in aa_orfs", aa_head

#   trans_orfs.each do |trans_head, trans_seq|
#     if aa_seq == trans_seq
#       tr_matches += 1
#       if ident_seqs.has_key? aa_head
#         ident_seqs[aa_head] << trans_head
#       else
#         uniq_tr_matches += 1
#         ident_seqs[aa_head] = [trans_head]
#       end
#     end
#   end
# end

Utils.log "Found %d translated matches (%d unique)",
          tr_matches,
          uniq_tr_matches

Utils.log "Found %d total matches", ident_seqs.count

still_need_nt_seqs = aa_orfs.keys - ident_seqs.keys

Utils.log "%d sequences did not have matches", still_need_nt_seqs.count

skip_count = 0
seq_bad = false
nt_aln_seq_lens = []

FastaFile.open(opts[:protein_alignment]).
  each_record do |aln_head, aln_seq|

  Utils.assert aa_orfs.has_key?(aln_head),
               "Header %s is present in file %s, but not in %s",
               aln_head,
               opts[:aa_orfs],
               opts[:protein_alignment]

  if ident_seqs.has_key? aln_head
    nt_headers = ident_seqs[aln_head]

    # TODO get to work for translated headers
    nt_headers.each do |nt_head|
      aln_seq_char_len = char_len(aln_seq)
      nt_orf_char_len = char_len(nt_orfs[nt_head])
      diff = aln_seq_char_len - nt_orf_char_len

      # TODO diff <= 3 cos not all AA have the stop at the end
      Utils.assert diff <= 3,
                   "Length mismatch for aln_seq: %s (%d) and " +
                   "nt_orf: %s (%d)",
                   aln_head,
                   aln_seq_char_len * 3,
                   nt_head,
                   nt_orf_char_len

      nt_orf = nt_orfs[nt_head].split ""
      new_nt_aln = ""
      aln_seq.each_char.with_index do |c, i|

        if i.zero? && c == "M" && nt_orf[0..2].join != "ATG"
          Utils.log "AA aln '%s' starts with 'M' but nt_orf '%s' " +
                    "starts with '%s'. " +
                    "AA has start, but NT doesn't.",
                    aln_head,
                    nt_head,
                    nt_orf[0..2].join
        elsif i.zero? && c != "M" && nt_orf[0..2].join == "ATG"
          Utils.log "AA aln '%s' starts with '%c' but nt_orf '%s' " +
                    "starts with '%s'. " +
                    "NT has start, but AA doesn't." ,
                    aln_head,
                    c,
                    nt_head,
                    nt_orf[0..2].join
        elsif AA_SET.include? c
          these_three = nt_orf.shift(3).join ""
          Utils.assert these_three.length == 3,
                       "Not enough bases in the NT orf '%s'",
                       nt_head

          unless AA_TO_DNA[c].include?(these_three)
            Utils.warn "AA aln '%s' (has '%c') and NT seq '%s' " +
                       "(has '%s') " +
                       "do not match at AA " +
                       "gapped pos %d. Skipping sequence.",
                       aln_head,
                       c,
                       nt_head,
                       these_three,
                       i + 1
            # TODO this isnt breaking at the right place
            skip_count += 1
            seq_bad = true
            break
          end

          new_nt_aln << these_three
        elsif c == "X"
          new_nt_aln << "N" * 3
        else
          # TODO should this just be c?
          new_nt_aln << "-" * 3
        end
      end

      unless seq_bad
        printf ">%s\n%s\n", aln_head, new_nt_aln
        nt_aln_seq_lens << new_nt_aln.length
      end
      seq_bad = false
    end
  end
end

Utils.log "%d sequences that had matches were skipped", skip_count

nt_aln_seq_len = nt_aln_seq_lens.first

Utils.log "Output alignment length: %d", nt_aln_seq_len
Utils.log "Num of seqs output: %d", nt_aln_seq_lens.count

Utils.assert nt_aln_seq_lens.all? { |len| len == nt_aln_seq_len },
             "DO NOT TRUST OUTPUT ALIGNMENT. " +
             "Output nucleotide alignments have different lengths"
