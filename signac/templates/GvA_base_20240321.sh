{% extends "base_script.sh" %}
{% block header %}
    {% block preamble %}
#!/bin/bash
#SBATCH --job-name="{{ id }}"
#SBATCH --account=def-hpcg1999
#SBATCH --qos=privileged
        {% set memory_requested = operations | calc_memory(parallel)  %}
        {% if memory_requested %}
#SBATCH --mem={{ memory_requested|format_memory }}
        {% endif %}
        {% if partition %}
#SBATCH --partition={{ partition }}
        {% endif %}
        {% set walltime = operations | calc_walltime(parallel) %}
        {% if walltime %}
#SBATCH -t {{ walltime|format_timedelta }}
        {% endif %}
        {% if job_output %}
#SBATCH --output={{ job_output }}
#SBATCH --error={{ job_output }}
        {% endif %}
    {% endblock preamble %}
{% endblock header %}
{% block custom_content %}
{#
    This block is not used by any other template and can be safely modified
    without the need to call super(). We recommend most additions to the
    templates go here if they are not direct changes to an existing template.

    For example, commands like `module load ...` or printing diagnostic
    information from the scheduler can be done in this block.
#}
module use /global/project/hpcg1999/software/modulefiles
module load GvA-base/2024.3.21
{% endblock custom_content %}
