{% extends "GvA_base_20240321.sh" %}

{% block header %}
    {{- super () -}}
#SBATCH -ntasks=20
#SBATCH -nodes=1
#SBATCH -cpus-per-task=1
#SBATCH --mem=2G
{% endblock header %}
{% block custom_content %}
{#
    This block is not used by any other template and can be safely modified
    without the need to call super(). We recommend most additions to the
    templates go here if they are not direct changes to an existing template.

    For example, commands like `module load ...` or printing diagnostic
    information from the scheduler can be done in this block.
#}
{% endblock custom_content %}
    {{- super () -}}
module load HypOptLib/0.1.1
{% block body %}
    {{- super () -}}
{% endblock body %}
