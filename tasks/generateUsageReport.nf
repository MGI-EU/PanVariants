process get_UsageReport{
    shell '/bin/bash', '-u'

    input:
        val usage_input

    script:
    """
    ${params.python3} ${params.generateUsageReport_py} ${params.outdir}/../
    cp ${params.outdir}/../runit_report.html ${params.outdir}/upload
    """
}
