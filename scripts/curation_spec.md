# Metadata Classification Script Specification (v1)

### Purpose

Generate one or more derived classification columns from a tab-delimited input file using a configurable rule file.

Typical use case:

```         
value   count
col-0   10792
columbia    3924
landsberg erecta    364
```

becomes

```         
value   count   col0    ler OR
col-0   10792   1   0   1
columbia    3924    1   0   1
landsberg erecta    364 0   1   1
```

## Input Data Format

The input file is tab-delimited.

Example:

```         
value   count
col-0   10792
columbia    3924
col0    752
landsberg erecta    364
```

### Requirements

-   First row is the header.
-   One column contains the text value to be classified.
-   Initial implementation assumes the first column is the classification target.
-   Rest columns are optional.

Example:

```         
value   count
```

where value is evaluated against the rules.

## Configuration File Format

The configuration file is tab-delimited.

### Columns

```         
column  type    rule
```

#### column

Output column name.

Examples:

```         
col0
ler
OR
```

#### type

Supported values:

```         
regex
logic
expr
```

#### rule

Definition of the rule.

Examples:

```         
col-0
^columbia$
col0 || ler
!OR
```

### Regex Rules

Regex rules accumulate matches into a column.

Example:

```         
column  type    rule
col0    regex   col-0
col0    regex   col_0
col0    regex   col0
col0    regex   ^columbia$
```

is logically equivalent to:

```         
$col0 ||= match('col-0');
$col0 ||= match('col_0');
$col0 ||= match('col0');
$col0 ||= match('^columbia$');
```

#### Behavior

A regex rule sets the target column to 1 if matched.

Once a column becomes 1, subsequent regex rules for the same column do not change it.

Effectively:

```         
$state{$column} ||= ($value =~ /$regex/i);
```

### Logic Rules

Logic rules perform direct assignment.

Example:

```         
OR  logic   col0 || ler
```

means:

```         
$OR = ($col0 || $ler) ? 1 : 0;
```

Example:

```         
NOT_OR  logic   !OR
```

means:

```         
$NOT_OR = (!$OR) ? 1 : 0;
```

#### Behavior

Logic rules assign the column value directly.

They can produce either:

```         
0
```

or

```         
1
```

The result of a logic rule replaces the previous value of that column.

Equivalent to:

```         
$state{$column} = evaluate_boolean_expression($rule);
```

### Expr Rules

Expr rules perform user-defined perl expressions. In the loop of processing each row, an array named `@fields` would be used to store the input fields. Expr rules are trusted Perl expressions evaluated with eval(). They may access `@fields` only. The result of the expression should be either 0 or 1, which will be assigned to the specified column.

Examples:

```
high_count	eval	$fields[1] > 1000

missing	eval	$fields[0] eq "missing"

common	eval	$fields[1] >= 100

interesting	eval	($fields[1] > 100) && ($fields[2] > 5)
```

### Evaluation Order

Rules are executed strictly in configuration file order.

Example:

```         
column  type    rule
col0    regex   col-0
ler regex   ler-0
OR  logic   col0 || ler
```

Execution:

```         
1. evaluate col0 regex
2. evaluate ler regex
3. evaluate OR logic
```

No separate regex phase or logic phase exists.

The configuration file defines the execution sequence.

### Variable Resolution

Logic expressions may reference previously defined columns.

Valid:

```         
col0    regex   col-0
ler regex   ler-0
OR  logic   col0 || ler
```

Invalid:

```         
OR  logic   col0 || ler
col0    regex   col-0
ler regex   ler-0
```

because col0 and ler have not yet been defined when OR is evaluated.

### Undefined Variables

Referencing an undefined column is a fatal error.

Example:

```         
OR  logic   col0 || ler
```

when ler has never been introduced earlier.

Error:

```         
ERROR: line 27
Undefined column 'ler' referenced in logic expression:
OR  logic   col0 || ler
```

Processing should terminate immediately.

### Supported Logic Syntax

Version 1 supports boolean expressions only.

#### Operators

```         
!
&&
||
(
)
```

Examples:

```         
col0
!col0

col0 && ler

col0 || ler

(col0 || ler) && WT

!(col0 || ler)
```

#### Constants

Supported:

```         
0
1
```

Examples:

```         
FOUND   logic   1

MISSING logic   0

MUT logic   !WT
```

### Regex Matching

Regex patterns are interpreted as Perl regular expressions.

Examples:

```         
^columbia$

efo_0005147

col-0
```

#### Case Sensitivity

Case-insensitive matching is handled by the script implementation.

The configuration file does not specify case sensitivity.

Current default behavior:

```         
$value =~ /$regex/i
```

### Output Columns

Output columns appear in the order of first appearance in the configuration file.

Example:

```         
col0    regex   col-0
ler regex   ler-0
OR  logic   col0 || ler
```

Output:

```         
value   count   col0    ler OR
```

Even if later rules reference the same columns again.

### Comments

Lines beginning with:

```         
#
```

are ignored.

Example:

```         
# Arabidopsis ecotypes

col0    regex   col-0
col0    regex   ^columbia$

# Landsberg erecta
ler regex   landsberg erecta

OR  logic   col0 || ler
```

## Example

Configuration:

```         
column  type    rule

col0    regex   col-0
col0    regex   col_0
col0    regex   ^columbia$

ler regex   landsberg erecta
ler regex   ler-0

OR  logic   col0 || ler
NOT_OR  logic   !OR
```

Input:

```         
value   count
columbia    100
landsberg erecta    200
unknown 300
```

Output:

```         
value   count   col0    ler OR  NOT_OR
columbia    100 1   0   1   0
landsberg erecta    200 0   1   1   0
unknown 300 0   0   0   1
```

This specification defines a simple rule engine where regex rules accumulate evidence and logic rules compute derived binary classifications, with execution strictly following the configuration file order.
