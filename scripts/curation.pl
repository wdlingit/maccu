#!/usr/bin/env perl

use strict;
use warnings;

my ($config_file, $input_file) = @ARGV;

die "Usage: $0 config.tsv input.tsv\n"
    unless defined $config_file && defined $input_file;

my (@operations, @output_columns);

read_config($config_file);

open(my $IN, "<", $input_file)
    or die "Cannot open $input_file: $!\n";

#
# Header
#
my $header = <$IN>;
defined $header or die "Empty input file\n";

chomp $header;

my @header = split(/\t/, $header);

print join("\t", @header, @output_columns), "\n";

#
# Data rows
#
while (<$IN>) {

    chomp;

    next if /^\s*$/;

    my @fields = split(/\t/, $_, -1);

    my $value = $fields[0];

    my %state;

    @state{@output_columns} = (0) x @output_columns;

    foreach my $op (@operations) {

        if ($op->{type} eq 'regex') {

            if ($value =~ $op->{regex}) {
                $state{$op->{column}} = 1;
            }

        }
        elsif ($op->{type} eq 'logic') {

            $state{$op->{column}}
                = $op->{code}->(\%state);
        }
        elsif ($op->{type} eq 'expr') {

            $state{$op->{column}}
                = $op->{code}->(\@fields);
        }
    }

    print join(
        "\t",
        @fields,
        @state{@output_columns}
    ), "\n";
}

close $IN;

##############################
# read config
##############################

sub read_config {

    my ($file) = @_;

    open(my $CFG, "<", $file)
        or die "Cannot open $file: $!\n";

    my %seen_column;
    my %known_column;

    my $line_no = 0;

    while (<$CFG>) {

        $line_no++;

        chomp;

        next if /^\s*$/;
        next if /^\s*#/;

        my ($column, $type, $rule) =
            split(/\t/, $_, 3);

        die "Config line $line_no: expected 3 columns\n"
            unless defined $column
            && defined $type
            && defined $rule;

        if (!$seen_column{$column}) {

            push @output_columns, $column;
            $seen_column{$column} = 1;
        }

        if ($type eq 'regex') {

            my $regex;

            eval {
                $regex = qr/$rule/i;
            };

            die "Config line $line_no: invalid regex '$rule'\n"
                if $@;

            push @operations, {
                line   => $line_no,
                column => $column,
                type   => 'regex',
                regex  => $regex,
            };
        }
        elsif ($type eq 'logic') {

            validate_logic(
                $rule,
                \%known_column,
                $line_no
            );

            my $code = compile_logic($rule);

            push @operations, {
                line   => $line_no,
                column => $column,
                type   => 'logic',
                code   => $code,
            };
        }
        elsif ($type eq 'expr') {
          
          validate_expr(
            $rule,
            $line_no
          );
          
          my $code =
            compile_expr($rule, $line_no);
            
          push @operations, {
            line   => $line_no,
            column => $column,
            type   => 'expr',
            code   => $code,
          };
        }
        else {

            die "Config line $line_no: unknown type '$type'\n";
        }

        #
        # Column becomes available after
        # its first appearance.
        #
        $known_column{$column} = 1;
    }

    close $CFG;
}

##############################
# validate expr expression
##############################

sub validate_expr {

    my ($expr, $line_no) = @_;

    #
    # forbid package access
    #
    if ($expr =~ /::/) {

        die
            "Config line $line_no: ".
            "package access (::) forbidden\n";
    }

    #
    # forbid shell execution
    #
    if ($expr =~ /`/) {

        die
            "Config line $line_no: ".
            "backticks forbidden\n";
    }

    #
    # obviously dangerous functions
    #
    if (
        $expr =~ /\b(
            system|
            exec|
            fork|
            open|
            unlink|
            rename|
            chdir|
            mkdir|
            rmdir|
            require|
            do|
            eval
        )\b/x
    ) {

        die
            "Config line $line_no: ".
            "forbidden function in expr\n";
    }

    #
    # only allow @fields references
    #
    while (
        $expr =~ /\$([A-Za-z_][A-Za-z0-9_]*)/g
    ) {

        my $var = $1;

        next if $var eq 'fields';

        die
            "Config line $line_no: ".
            "expr may only reference ".
            "\$fields[] (found \$$var)\n";
    }
}

##############################
# compile expr
##############################

sub compile_expr {

    my ($expr, $line_no) = @_;

    my $perl_expr = $expr;

    #
    # Convert:
    #   $fields[1]
    # into:
    #   $fields->[1]
    #
    $perl_expr =~ s/
        \$fields
        \[
        (\d+)
        \]
    /
        \$fields->[$1]
    /gx;

    my $source = qq{
        sub {
            my (\$fields) = \@_;

            return (
                ($perl_expr)
                ? 1
                : 0
            );
        }
    };

    my $sub = eval $source;

    if ($@) {

        die
            "Config line $line_no: ".
            "cannot compile expr:\n".
            "$expr\n$@";
    }

    return $sub;
}

##############################
# validate logic expression
##############################

sub validate_logic {

    my ($expr, $known, $line_no) = @_;

    while ($expr =~ /\b([A-Za-z_][A-Za-z0-9_]*)\b/g) {

        my $token = $1;

        die join(
            "",
            "Config line $line_no: ",
            "undefined column '$token'\n"
        )
        unless exists $known->{$token};
    }
}

##############################
# compile logic
##############################

sub compile_logic {

    my ($expr) = @_;

    my $perl_expr = $expr;

    $perl_expr =~ s{
        \b([A-Za-z_][A-Za-z0-9_]*)\b
    }{
        "\$state->{'$1'}"
    }gex;

    my $source = qq{
        sub {
            my \$state = shift;
            return (($perl_expr) ? 1 : 0);
        }
    };

    my $sub = eval $source;

    die "Cannot compile logic expression:\n$expr\n$@\n"
        if $@;

    return $sub;
}
