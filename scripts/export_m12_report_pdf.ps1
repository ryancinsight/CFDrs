param(
    [string]$InputPath = ".\\report\\ARPA-H_SonALAsense_Milestone 12 Report.md",
    [string]$OutputPath = ".\\report\\ARPA-H_SonALAsense_Milestone 12 Report.pdf"
)

$ErrorActionPreference = "Stop"

function Convert-LatexFragmentToText {
    param([string]$Text)

    $converted = $Text
    $converted = $converted -replace '\\frac\s*\{([^{}]+)\}\s*\{([^{}]+)\}', '($1)/($2)'
    $converted = $converted -replace '\\text\s*\{([^{}]+)\}', '$1'
    $converted = $converted -replace '\\mathrm\s*\{([^{}]+)\}', '$1'
    $converted = $converted -replace '\\mathbf\s*\{([^{}]+)\}', '$1'
    $converted = $converted -replace '\\operatorname\s*\{([^{}]+)\}', '$1'

    $replacements = @{
        '\\propto' = ' proportional to '
        '\\approx' = ' approx '
        '\\times' = ' x '
        '\\cdot' = ' * '
        '\\mu' = 'u'
        '\\kappa' = 'kappa'
        '\\eta' = 'eta'
        '\\sigma' = 'sigma'
        '\\Sigma' = 'Sigma'
        '\\sum' = 'sum'
        '\\Delta' = 'Delta'
        '\\delta' = 'delta'
        '\\alpha' = 'alpha'
        '\\beta' = 'beta'
        '\\gamma' = 'gamma'
        '\\pi' = 'pi'
        '\\lambda' = 'lambda'
        '\\infty' = 'infinity'
        '\\leq' = '<='
        '\\geq' = '>='
        '\\neq' = '!='
        '\\to' = '->'
        '\\Rightarrow' = '=>'
        '\\left' = ''
        '\\right' = ''
    }

    foreach ($key in $replacements.Keys) {
        $converted = $converted -replace $key, $replacements[$key]
    }

    $converted = $converted -replace '\\_', '_'
    $converted = $converted -replace '\\\^', '^'
    $converted = $converted -replace '\\([{}])', '$1'
    $converted = $converted -replace '\s+', ' '
    return $converted.Trim()
}

function Convert-MarkdownForPdf {
    param([string]$Markdown)

        $tableStyle = @"
<style>
table {
    width: 100% !important;
    table-layout: fixed;
    font-size: 8.5px;
}

th, td {
    font-size: 8.5px;
    padding: 3px 5px;
    vertical-align: top;
    white-space: normal;
    overflow-wrap: anywhere;
    word-break: break-word;
}

table code {
    font-size: 7.5px;
    white-space: pre-wrap;
    overflow-wrap: anywhere;
    word-break: break-word;
}
</style>

"@

        $converted = $tableStyle + ($Markdown -replace '</?mark[^>]*>', '')

    $blockPattern = '(?s)\$\$(.+?)\$\$'
    $converted = [regex]::Replace($converted, $blockPattern, {
        param($match)
        $body = Convert-LatexFragmentToText $match.Groups[1].Value
        "`r`n```text`r`n$body`r`n```"
    })

    $inlinePattern = '(?<!\$)\$([^`\r\n]+?)\$(?!\$)'
    $converted = [regex]::Replace($converted, $inlinePattern, {
        param($match)
        Convert-LatexFragmentToText $match.Groups[1].Value
    })

    return $converted
}

$resolvedInput = Resolve-Path $InputPath
$resolvedOutput = [System.IO.Path]::GetFullPath($OutputPath)
$tempPath = Join-Path ([System.IO.Path]::GetTempPath()) ([System.IO.Path]::GetRandomFileName() + '.md')

try {
    $markdown = Get-Content -LiteralPath $resolvedInput -Raw
    $pdfSafeMarkdown = Convert-MarkdownForPdf $markdown
    Set-Content -LiteralPath $tempPath -Value $pdfSafeMarkdown -Encoding UTF8

    md-to-pdf $tempPath --dest $resolvedOutput
}
finally {
    if (Test-Path -LiteralPath $tempPath) {
        Remove-Item -LiteralPath $tempPath -Force
    }
}