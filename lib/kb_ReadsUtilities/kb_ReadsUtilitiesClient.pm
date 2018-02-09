package kb_ReadsUtilities::kb_ReadsUtilitiesClient;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

kb_ReadsUtilities::kb_ReadsUtilitiesClient

=head1 DESCRIPTION


** A KBase module: kb_ReadsUtilities
**
** This module contains basic utility Apps for manipulating Reads Libraries


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => kb_ReadsUtilities::kb_ReadsUtilitiesClient::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my %arg_hash2 = @args;
	if (exists $arg_hash2{"token"}) {
	    $self->{token} = $arg_hash2{"token"};
	} elsif (exists $arg_hash2{"user_id"}) {
	    my $token = Bio::KBase::AuthToken->new(@args);
	    if (!$token->error_message) {
	        $self->{token} = $token->token;
	    }
	}
	
	if (exists $self->{token})
	{
	    $self->{client}->{token} = $self->{token};
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 KButil_FASTQ_to_FASTA

  $return = $obj->KButil_FASTQ_to_FASTA($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_ReadsUtilities.KButil_FASTQ_to_FASTA_Params
$return is a kb_ReadsUtilities.KButil_FASTQ_to_FASTA_Output
KButil_FASTQ_to_FASTA_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
	output_name has a value which is a kb_ReadsUtilities.data_obj_name
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_FASTQ_to_FASTA_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_ReadsUtilities.KButil_FASTQ_to_FASTA_Params
$return is a kb_ReadsUtilities.KButil_FASTQ_to_FASTA_Output
KButil_FASTQ_to_FASTA_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
	output_name has a value which is a kb_ReadsUtilities.data_obj_name
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_FASTQ_to_FASTA_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_FASTQ_to_FASTA
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_FASTQ_to_FASTA (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_FASTQ_to_FASTA:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_FASTQ_to_FASTA');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_ReadsUtilities.KButil_FASTQ_to_FASTA",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_FASTQ_to_FASTA',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_FASTQ_to_FASTA",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_FASTQ_to_FASTA',
				       );
    }
}
 


=head2 KButil_Split_Reads

  $return = $obj->KButil_Split_Reads($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_ReadsUtilities.KButil_Split_Reads_Params
$return is a kb_ReadsUtilities.KButil_Split_Reads_Output
KButil_Split_Reads_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
	output_name has a value which is a kb_ReadsUtilities.data_obj_name
	split_num has a value which is an int
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Split_Reads_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_ReadsUtilities.KButil_Split_Reads_Params
$return is a kb_ReadsUtilities.KButil_Split_Reads_Output
KButil_Split_Reads_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
	output_name has a value which is a kb_ReadsUtilities.data_obj_name
	split_num has a value which is an int
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Split_Reads_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Split_Reads
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Split_Reads (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Split_Reads:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Split_Reads');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_ReadsUtilities.KButil_Split_Reads",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Split_Reads',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Split_Reads",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Split_Reads',
				       );
    }
}
 


=head2 KButil_Random_Subsample_Reads

  $return = $obj->KButil_Random_Subsample_Reads($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_ReadsUtilities.KButil_Random_Subsample_Reads_Params
$return is a kb_ReadsUtilities.KButil_Random_Subsample_Reads_Output
KButil_Random_Subsample_Reads_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
	output_name has a value which is a kb_ReadsUtilities.data_obj_name
	subsample_fraction has a value which is a kb_ReadsUtilities.Fractionate_Options
	desc has a value which is a string
	seed has a value which is an int
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
Fractionate_Options is a reference to a hash where the following keys are defined:
	split_num has a value which is an int
	reads_num has a value which is an int
	reads_perc has a value which is a float
KButil_Random_Subsample_Reads_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_ReadsUtilities.KButil_Random_Subsample_Reads_Params
$return is a kb_ReadsUtilities.KButil_Random_Subsample_Reads_Output
KButil_Random_Subsample_Reads_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
	output_name has a value which is a kb_ReadsUtilities.data_obj_name
	subsample_fraction has a value which is a kb_ReadsUtilities.Fractionate_Options
	desc has a value which is a string
	seed has a value which is an int
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
Fractionate_Options is a reference to a hash where the following keys are defined:
	split_num has a value which is an int
	reads_num has a value which is an int
	reads_perc has a value which is a float
KButil_Random_Subsample_Reads_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Random_Subsample_Reads
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Random_Subsample_Reads (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Random_Subsample_Reads:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Random_Subsample_Reads');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_ReadsUtilities.KButil_Random_Subsample_Reads",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Random_Subsample_Reads',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Random_Subsample_Reads",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Random_Subsample_Reads',
				       );
    }
}
 


=head2 KButil_Merge_ReadsSet_to_OneLibrary

  $return = $obj->KButil_Merge_ReadsSet_to_OneLibrary($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_ReadsUtilities.KButil_Merge_ReadsSet_to_OneLibrary_Params
$return is a kb_ReadsUtilities.KButil_Merge_ReadsSet_to_OneLibrary_Output
KButil_Merge_ReadsSet_to_OneLibrary_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
	output_name has a value which is a kb_ReadsUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Merge_ReadsSet_to_OneLibrary_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_ReadsUtilities.KButil_Merge_ReadsSet_to_OneLibrary_Params
$return is a kb_ReadsUtilities.KButil_Merge_ReadsSet_to_OneLibrary_Output
KButil_Merge_ReadsSet_to_OneLibrary_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
	output_name has a value which is a kb_ReadsUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Merge_ReadsSet_to_OneLibrary_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Merge_ReadsSet_to_OneLibrary
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Merge_ReadsSet_to_OneLibrary (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Merge_ReadsSet_to_OneLibrary:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Merge_ReadsSet_to_OneLibrary');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_ReadsUtilities.KButil_Merge_ReadsSet_to_OneLibrary",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Merge_ReadsSet_to_OneLibrary',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Merge_ReadsSet_to_OneLibrary",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Merge_ReadsSet_to_OneLibrary',
				       );
    }
}
 


=head2 KButil_Merge_MultipleReadsLibs_to_OneLibrary

  $return = $obj->KButil_Merge_MultipleReadsLibs_to_OneLibrary($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_ReadsUtilities.KButil_Merge_MultipleReadsLibs_to_OneLibrary_Params
$return is a kb_ReadsUtilities.KButil_Merge_MultipleReadsLibs_to_OneLibrary_Output
KButil_Merge_MultipleReadsLibs_to_OneLibrary_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_refs has a value which is a kb_ReadsUtilities.data_obj_ref
	output_name has a value which is a kb_ReadsUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Merge_MultipleReadsLibs_to_OneLibrary_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_ReadsUtilities.KButil_Merge_MultipleReadsLibs_to_OneLibrary_Params
$return is a kb_ReadsUtilities.KButil_Merge_MultipleReadsLibs_to_OneLibrary_Output
KButil_Merge_MultipleReadsLibs_to_OneLibrary_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_refs has a value which is a kb_ReadsUtilities.data_obj_ref
	output_name has a value which is a kb_ReadsUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Merge_MultipleReadsLibs_to_OneLibrary_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Merge_MultipleReadsLibs_to_OneLibrary
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Merge_MultipleReadsLibs_to_OneLibrary (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Merge_MultipleReadsLibs_to_OneLibrary:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Merge_MultipleReadsLibs_to_OneLibrary');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_ReadsUtilities.KButil_Merge_MultipleReadsLibs_to_OneLibrary",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Merge_MultipleReadsLibs_to_OneLibrary',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Merge_MultipleReadsLibs_to_OneLibrary",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Merge_MultipleReadsLibs_to_OneLibrary',
				       );
    }
}
 


=head2 KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs

  $return = $obj->KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_ReadsUtilities.KButil_Extract_Unpaired_Reads_Params
$return is a kb_ReadsUtilities.KButil_Extract_Unpaired_Reads_Output
KButil_Extract_Unpaired_Reads_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
	output_name has a value which is a kb_ReadsUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Extract_Unpaired_Reads_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_ReadsUtilities.KButil_Extract_Unpaired_Reads_Params
$return is a kb_ReadsUtilities.KButil_Extract_Unpaired_Reads_Output
KButil_Extract_Unpaired_Reads_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
	output_name has a value which is a kb_ReadsUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Extract_Unpaired_Reads_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_ReadsUtilities.KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs',
				       );
    }
}
 


=head2 KButil_Translate_ReadsLibs_QualScores

  $return = $obj->KButil_Translate_ReadsLibs_QualScores($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_ReadsUtilities.KButil_Translate_ReadsLibs_QualScores_Params
$return is a kb_ReadsUtilities.KButil_Translate_ReadsLibs_QualScores_Output
KButil_Translate_ReadsLibs_QualScores_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_refs has a value which is a kb_ReadsUtilities.data_obj_ref
workspace_name is a string
data_obj_ref is a string
KButil_Translate_ReadsLibs_QualScores_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref
data_obj_name is a string

</pre>

=end html

=begin text

$params is a kb_ReadsUtilities.KButil_Translate_ReadsLibs_QualScores_Params
$return is a kb_ReadsUtilities.KButil_Translate_ReadsLibs_QualScores_Output
KButil_Translate_ReadsLibs_QualScores_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_refs has a value which is a kb_ReadsUtilities.data_obj_ref
workspace_name is a string
data_obj_ref is a string
KButil_Translate_ReadsLibs_QualScores_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref
data_obj_name is a string


=end text

=item Description



=back

=cut

 sub KButil_Translate_ReadsLibs_QualScores
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Translate_ReadsLibs_QualScores (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Translate_ReadsLibs_QualScores:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Translate_ReadsLibs_QualScores');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_ReadsUtilities.KButil_Translate_ReadsLibs_QualScores",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Translate_ReadsLibs_QualScores',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Translate_ReadsLibs_QualScores",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Translate_ReadsLibs_QualScores',
				       );
    }
}
 


=head2 KButil_AddInsertLen_to_ReadsLibs

  $return = $obj->KButil_AddInsertLen_to_ReadsLibs($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_ReadsUtilities.KButil_AddInsertLen_to_ReadsLibs_Params
$return is a kb_ReadsUtilities.KButil_AddInsertLen_to_ReadsLibs_Output
KButil_AddInsertLen_to_ReadsLibs_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_refs has a value which is a kb_ReadsUtilities.data_obj_ref
	insert_len has a value which is a float
	insert_stddev has a value which is a float
workspace_name is a string
data_obj_ref is a string
KButil_AddInsertLen_to_ReadsLibs_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref
data_obj_name is a string

</pre>

=end html

=begin text

$params is a kb_ReadsUtilities.KButil_AddInsertLen_to_ReadsLibs_Params
$return is a kb_ReadsUtilities.KButil_AddInsertLen_to_ReadsLibs_Output
KButil_AddInsertLen_to_ReadsLibs_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_ReadsUtilities.workspace_name
	input_refs has a value which is a kb_ReadsUtilities.data_obj_ref
	insert_len has a value which is a float
	insert_stddev has a value which is a float
workspace_name is a string
data_obj_ref is a string
KButil_AddInsertLen_to_ReadsLibs_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_ReadsUtilities.data_obj_name
	report_ref has a value which is a kb_ReadsUtilities.data_obj_ref
data_obj_name is a string


=end text

=item Description



=back

=cut

 sub KButil_AddInsertLen_to_ReadsLibs
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_AddInsertLen_to_ReadsLibs (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_AddInsertLen_to_ReadsLibs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_AddInsertLen_to_ReadsLibs');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_ReadsUtilities.KButil_AddInsertLen_to_ReadsLibs",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_AddInsertLen_to_ReadsLibs',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_AddInsertLen_to_ReadsLibs",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_AddInsertLen_to_ReadsLibs',
				       );
    }
}
 
  
sub status
{
    my($self, @args) = @_;
    if ((my $n = @args) != 0) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function status (received $n, expecting 0)");
    }
    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
        method => "kb_ReadsUtilities.status",
        params => \@args,
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => 'status',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
                          );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method status",
                        status_line => $self->{client}->status_line,
                        method_name => 'status',
                       );
    }
}
   

sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_ReadsUtilities.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'KButil_AddInsertLen_to_ReadsLibs',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method KButil_AddInsertLen_to_ReadsLibs",
            status_line => $self->{client}->status_line,
            method_name => 'KButil_AddInsertLen_to_ReadsLibs',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for kb_ReadsUtilities::kb_ReadsUtilitiesClient\n";
    }
    if ($sMajor == 0) {
        warn "kb_ReadsUtilities::kb_ReadsUtilitiesClient version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 workspace_name

=over 4



=item Description

** The workspace object refs are of form:
**
**    objects = ws.get_objects([{'ref': params['workspace_id']+'/'+params['obj_name']}])
**
** "ref" means the entire name combining the workspace id and the object name
** "id" is a numerical identifier of the workspace or object, and should just be used for workspace
** "name" is a string identifier of a workspace or object.  This is received from Narrative.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 sequence

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 data_obj_name

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 data_obj_ref

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 bool

=over 4



=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 KButil_FASTQ_to_FASTA_Params

=over 4



=item Description

KButil_FASTQ_to_FASTA()
**
** Method for Converting a FASTQ SingleEndLibrary to a FASTA SingleEndLibrary


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
output_name has a value which is a kb_ReadsUtilities.data_obj_name

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
output_name has a value which is a kb_ReadsUtilities.data_obj_name


=end text

=back



=head2 KButil_FASTQ_to_FASTA_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=back



=head2 KButil_Split_Reads_Params

=over 4



=item Description

KButil_Split_Reads()
**
**  Method for spliting a ReadsLibrary into evenly sized ReadsLibraries


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
output_name has a value which is a kb_ReadsUtilities.data_obj_name
split_num has a value which is an int
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
output_name has a value which is a kb_ReadsUtilities.data_obj_name
split_num has a value which is an int
desc has a value which is a string


=end text

=back



=head2 KButil_Split_Reads_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=back



=head2 Fractionate_Options

=over 4



=item Description

KButil_Random_Subsample_Reads()
**
**  Method for random subsampling of reads library


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
split_num has a value which is an int
reads_num has a value which is an int
reads_perc has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
split_num has a value which is an int
reads_num has a value which is an int
reads_perc has a value which is a float


=end text

=back



=head2 KButil_Random_Subsample_Reads_Params

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
output_name has a value which is a kb_ReadsUtilities.data_obj_name
subsample_fraction has a value which is a kb_ReadsUtilities.Fractionate_Options
desc has a value which is a string
seed has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
output_name has a value which is a kb_ReadsUtilities.data_obj_name
subsample_fraction has a value which is a kb_ReadsUtilities.Fractionate_Options
desc has a value which is a string
seed has a value which is an int


=end text

=back



=head2 KButil_Random_Subsample_Reads_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=back



=head2 KButil_Merge_ReadsSet_to_OneLibrary_Params

=over 4



=item Description

KButil_Merge_ReadsSet_to_OneLibrary()
**
**  Method for merging a ReadsSet into one library


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
output_name has a value which is a kb_ReadsUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
output_name has a value which is a kb_ReadsUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Merge_ReadsSet_to_OneLibrary_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=back



=head2 KButil_Merge_MultipleReadsLibs_to_OneLibrary_Params

=over 4



=item Description

KButil_Merge_MultipleReadsLibs_to_OneLibrary()
**
**  Method for merging ReadsLibs into one library


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_refs has a value which is a kb_ReadsUtilities.data_obj_ref
output_name has a value which is a kb_ReadsUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_refs has a value which is a kb_ReadsUtilities.data_obj_ref
output_name has a value which is a kb_ReadsUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Merge_MultipleReadsLibs_to_OneLibrary_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=back



=head2 KButil_Extract_Unpaired_Reads_Params

=over 4



=item Description

KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs()
**
**  Method for removing unpaired reads from a paired end library or set and matching the order of reads


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
output_name has a value which is a kb_ReadsUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_ref has a value which is a kb_ReadsUtilities.data_obj_ref
output_name has a value which is a kb_ReadsUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Extract_Unpaired_Reads_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=back



=head2 KButil_Translate_ReadsLibs_QualScores_Params

=over 4



=item Description

KButil_Translate_ReadsLibs_QualScores()
**
**  Method for Translating ReadsLibs Qual scores


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_refs has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_refs has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=back



=head2 KButil_Translate_ReadsLibs_QualScores_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=back



=head2 KButil_AddInsertLen_to_ReadsLibs_Params

=over 4



=item Description

KButil_AddInsertLen_to_ReadsLibs()
**
**  Method for Adding Insert Len to PairedEnd ReadsLibs


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_refs has a value which is a kb_ReadsUtilities.data_obj_ref
insert_len has a value which is a float
insert_stddev has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_ReadsUtilities.workspace_name
input_refs has a value which is a kb_ReadsUtilities.data_obj_ref
insert_len has a value which is a float
insert_stddev has a value which is a float


=end text

=back



=head2 KButil_AddInsertLen_to_ReadsLibs_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_ReadsUtilities.data_obj_name
report_ref has a value which is a kb_ReadsUtilities.data_obj_ref


=end text

=back



=cut

package kb_ReadsUtilities::kb_ReadsUtilitiesClient::RpcClient;
use base 'JSON::RPC::Client';
use POSIX;
use strict;

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
