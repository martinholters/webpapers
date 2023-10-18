const make_example_table_data = (data) => {
  const f = (x) => -Math.pow(x, 3);
  const df = (x) => -3 * Math.pow(x, 2);
  const result = Array.from({ length: 6 }, (_, n) => {
    return { n };
  });
  // TR
  let x = 1;
  for (let n = 0; n < 6; n++) {
    result[n].trapz = x;
    let xn = x;
    for (let i = 0; i < 10; i++) {
      let r = x + 0.5 * (f(xn) + f(x)) - xn;
      xn = xn - r / (0.5 * df(xn) - 1);
    }
    x = xn;
  }
  // TR alternative
  x = 1 + 0.5 * f(1);
  for (let n = 0; n < 6; n++) {
    result[n].trapzalt = x;
    let xn = x;
    for (let i = 0; i < 10; i++) {
      let r = x + f(0.5 * (xn + x)) - xn;
      xn = xn - r / (0.5 * df(0.5 * (xn + x)) - 1);
    }
    x = xn;
  }
  // IMP
  x = 1;
  for (let n = 0; n < 6; n++) {
    result[n].imp = x;
    let xn = x;
    for (let i = 0; i < 10; i++) {
      let r = x + f(0.5 * (xn + x)) - xn;
      xn = xn - r / (0.5 * df(0.5 * (xn + x)) - 1);
    }
    x = xn;
  }
  // IMP alternative
  x = 0.5 * (1 + result[1].imp);
  for (let n = 0; n < 6; n++) {
    result[n].impalt = x;
    let xn = x;
    for (let i = 0; i < 10; i++) {
      let r = x + 0.5 * (f(xn) + f(x)) - xn;
      xn = xn - r / (0.5 * df(xn) - 1);
    }
    x = xn;
  }
  return result;
};

module.exports = {
  eleventyComputed: {
    example_results: make_example_table_data,
  },
};
