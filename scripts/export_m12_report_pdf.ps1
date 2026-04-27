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

    $replacements = @(
        @('\\propto', ' proportional to '),
        @('\\approx', ' approx '),
        @('\\times', ' x '),
        @('\\cdot', ' * '),
        @('\\mu', 'u'),
        @('\\kappa', 'kappa'),
        @('\\eta', 'eta'),
        @('\\sigma', 'sigma'),
        @('\\Sigma', 'Sigma'),
        @('\\sum', 'sum'),
        @('\\Delta', 'Delta'),
        @('\\delta', 'delta'),
        @('\\alpha', 'alpha'),
        @('\\beta', 'beta'),
        @('\\gamma', 'gamma'),
        @('\\pi', 'pi'),
        @('\\lambda', 'lambda'),
        @('\\infty', 'infinity'),
        @('\\leq', '<='),
        @('\\geq', '>='),
        @('\\neq', '!='),
        @('\\to', '->'),
        @('\\Rightarrow', '=>'),
        @('\\left', ''),
        @('\\right', '')
    )

    foreach ($replacement in $replacements) {
        $converted = $converted -replace $replacement[0], $replacement[1]
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

    $unicodeReplacements = @(
        @([string][char]0x03A3, 'sum'),
        @([string][char]0x1D62, '_i'),
        @([string][char]0x0394, 'Delta'),
        @([string][char]0x03B4, 'delta'),
        @([string][char]0x03C3, 'sigma'),
        @([string][char]0x03BA, 'kappa'),
        @([string][char]0x03B7, 'eta'),
        @([string][char]0x03BC, 'u'),
        @([string][char]0x00B5, 'u'),
        @([string][char]0x03C0, 'pi'),
        @([string][char]0x2192, '->'),
        @([string][char]0x2264, '<='),
        @([string][char]0x2265, '>='),
        @([string][char]0x00B1, '+/-')
    )
    foreach ($replacement in $unicodeReplacements) {
        $converted = $converted -replace $replacement[0], $replacement[1]
    }

    $blockPattern = '(?s)\$\$(.+?)\$\$'
    $converted = [regex]::Replace($converted, $blockPattern, {
        param($match)
        $body = Convert-LatexFragmentToText $match.Groups[1].Value
        "`r`n``````text`r`n$body`r`n``````"
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
$inputDirectory = Split-Path -Parent $resolvedInput
$tempPath = Join-Path $inputDirectory ('.m12_pdf_' + [System.Guid]::NewGuid().ToString('N') + '.md')
$tempPdfPath = [System.IO.Path]::ChangeExtension($tempPath, '.pdf')

try {
    $markdown = Get-Content -LiteralPath $resolvedInput -Raw -Encoding UTF8
    $pdfSafeMarkdown = Convert-MarkdownForPdf $markdown
    Set-Content -LiteralPath $tempPath -Value $pdfSafeMarkdown -Encoding UTF8

    npx --yes --package md-to-pdf -- md-to-pdf $tempPath
    if ($LASTEXITCODE -ne 0) {
        throw "md-to-pdf failed with exit code $LASTEXITCODE"
    }
    if (-not (Test-Path -LiteralPath $tempPdfPath)) {
        throw "md-to-pdf did not create the expected PDF at $tempPdfPath"
    }
    if ((Get-Item -LiteralPath $tempPdfPath).Length -le 0) {
        throw "md-to-pdf created an empty PDF at $tempPdfPath"
    }
    Move-Item -LiteralPath $tempPdfPath -Destination $resolvedOutput -Force
}
finally {
    if (Test-Path -LiteralPath $tempPath) {
        Remove-Item -LiteralPath $tempPath -Force
    }
    if (Test-Path -LiteralPath $tempPdfPath) {
        Remove-Item -LiteralPath $tempPdfPath -Force
    }
}
