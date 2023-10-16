const { EleventyHtmlBasePlugin } = require("@11ty/eleventy");
const katex = require('katex');
const cheerio = require('cheerio');

const mathrenderer = (md) => {
  const orig_fence_renderer = md.renderer.rules.fence;
  md.renderer.rules.fence = (...args) => {
    const [tokens, n] = args;
    const page = args[3].page;
    if (!page.mathinfo) {
      page.mathinfo = { eqctr: 1, labels: {} };
    }
    const isEq = /^equation\s?/.test(tokens[n].info);
    if (!isEq && tokens[n].info != 'math') {
      return orig_fence_renderer(...args);
    }
    const eqTag = isEq ? `\\tag{${page.mathinfo.eqctr}}` : '';
    const label = tokens[n].info.match(/equation\s+label=(.*)/)?.[1];
    if (label) {
      page.mathinfo.labels[label] = page.mathinfo.eqctr;
    }
    const html = katex.renderToString(`${eqTag}${tokens[n].content}`, {
      displayMode: true,
    });
    if (isEq) {
      page.mathinfo.eqctr++;
      if (label) {
        return `<span id="eq:${label}">${html}</span>`;
      }
    }
    return html;
  };
  const orig_code_inline_renderer = md.renderer.rules.code_inline;
  md.renderer.rules.code_inline = (...args) => {
    const [tokens, n] = args;
    const content = tokens[n].content;
    if (content[0] != '$' || content[content.length - 1] != '$') {
      return orig_code_inline_renderer(...args);
    }
    return katex.renderToString(content.substring(1, content.length - 1), {
      displayMode: false,
    });
  };
};

module.exports = function (eleventyConfig) {
  eleventyConfig.addPlugin(EleventyHtmlBasePlugin);
  eleventyConfig.addPassthroughCopy({ 'node_modules/katex/dist/fonts/*': 'fonts' });
  eleventyConfig.amendLibrary('md', (md) => {
    md.use(mathrenderer);
  });
  eleventyConfig.addShortcode('eqref', function (label) {
    return `(<a data-eqref="${label}">[unresolved reference]</a>)`;
  });
  eleventyConfig.addFilter('fixrefs', function (value) {
    const dom = cheerio.load(value, null, false);
    for (const el of Array.from(dom('a[data-eqref]'), dom)) {
      const label = el.attr('data-eqref');
      el.attr('href', `#eq:${label}`);
      el.text(`${this.page.mathinfo.labels[label]}`);
    }
    return dom.html();
  });
};
